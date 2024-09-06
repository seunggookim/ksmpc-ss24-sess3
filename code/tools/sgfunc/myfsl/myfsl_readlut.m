function tblout = myfsl_readlut(fnlut)
% tblout = myfsl_readlut(fnlut)
%
%
% [fslpath,'/fslpython/envs/fslpython/lib/python*/',...
%   'site-packages/fsleyes/assets/luts/harvard-oxford-cortical.lut']
% based on the fixed columns of strings.

tbl = readtable(fnlut,'FileType','text', 'delimiter','\t',...
  'NumHeaderLines',0, 'ReadVariableNames',false);
tblout = [];
for i = 1:numel(tbl)
  tblout.label(i,1) = str2double(tbl.Var1{i}(1:2));
  tblout.r(i,1) = str2double(tbl.Var1{i}(4:10));
  tblout.g(i,1) = str2double(tbl.Var1{i}(12:18));
  tblout.b(i,1) = str2double(tbl.Var1{i}(20:26));
  tblout.name{i,1} = tbl.Var1{i}(28:end);
end
tblout = struct2table(tblout);

end