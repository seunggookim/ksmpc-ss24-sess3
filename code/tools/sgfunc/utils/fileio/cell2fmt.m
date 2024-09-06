function fmtstr = cell2fmt (fmtcell, delimiter)
% fmtstr = cell2fmt (fmtcell, delimiter)
%
% (cc) sgKIM
fmtstr='';
if ~exist('delimiter','var')
delimiter='\t';
end

for i=1:numel(fmtcell)
fmtstr=[fmtstr delimiter fmtcell{i}];
end
fmtstr(1:2)='';
fmtstr=[fmtstr '\n'];
end
