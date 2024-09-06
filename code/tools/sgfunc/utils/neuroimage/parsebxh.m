function bxh = parsebxh(fn_bxh)
% bxh = parsebxh(fn_bxh)
%
% parse XML BIAC header
% (cc) 2019, sgKIM

xml = xmlread(fn_bxh);
elemnode = xml.getElementsByTagName('acquisitiondata');
acqdata = elemnode.item(0);
k = acqdata.getLength;
bxh = [];
for i=0:k-1
  fieldname = char(acqdata.item(i).getNodeName);
  value = char(acqdata.item(i).getTextContent);
  if ~strcmp(fieldname(1),'#')
    if ~isnan(str2num(value)) % str2num is used for a reason!
      value = str2num(value);
    end
    bxh.(fieldname) = value;
  end
end

% slice timing order
elemnode = xml.getElementsByTagName('datapoints');
textnode = elemnode.item(0); %  Currently the first datapoints tag is for
% label="acquisitiontimeindex"
bxh.acquisitiontimeindex = str2num(char(textnode.getTextContent));

% slice timing in msec, assuming # slices per band was an odd-multiple of
% the # of bands...
slice_per_band = max(bxh.acquisitiontimeindex);
% sanity check:
if slice_per_band * bxh.mb_factor ~= numel(bxh.acquisitiontimeindex)
  warning('# of slices per band was NOT a multiple of the # of bands!')
end
bxh.acquisitiontime_msec = (bxh.acquisitiontimeindex-1)*(bxh.tr/slice_per_band);

end