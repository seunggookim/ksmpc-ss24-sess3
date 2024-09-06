function [p1,f1,e1]=myfileparts(fname)
% [p1,f1,e1]=myfileparts(fname)
% (cc) sgKIM.
sp = filesep;
[p1,f1,e1] = fileparts(fname);
if isempty(p1) % if not pathname is given
 p1=pwd;
end
if strcmp(p1(1),'~') % if it was the home directory
  if ispc
    p1 = [getenv('USERPROFILE') p1(2:end)];
    % REF: https://github.com/geospace-code/matlab-stdlib/blob/main/%2Bstdlib/%2Bfileio/expanduser.m
  else
    p1 = [getenv('HOME') p1(2:end)];
  end

end
if ~strcmp(p1(1),sp) % if the path is relative
 p1=[pwd,sp,p1];
end
if strcmp(e1,'.gz') % if the filename is nii.gz or gii.gz or similar
 e1=[f1(end-3:end),e1];
 f1=f1(1:end-4);
end

end