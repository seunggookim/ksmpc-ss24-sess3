function impedance = myeeg_bvimp(fn_vhdr)
% import impedance readings from BrainVision Header (VHDR) file
%
% impedance = myeeg_bvimp(fn_vhdr)
%
% (cc) 2021, sgKIM.

fid = fopen(fn_vhdr,'r');
thisline = '';
while ~contains(thisline, 'Impedance [KOhm]')
  thisline = fgetl(fid);
end
C = textscan(fid,'%d%s%s','Delimiter',':');
C{3} = cellfun(@str2num, C{3}, 'uni', false);
ind = ~cellfun(@numel,C{3});
if sum(ind)
  [C{3}{ind}] = deal(inf); % "Out of Range!"
end
impedance = struct(...
  'chanindex',num2cell(C{1}), 'channame',C{2}, 'KOhm',C{3});
end