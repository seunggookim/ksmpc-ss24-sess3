function fn_nii = myspm_read_gz(fn)
% fn_nii = myspm_read_gz(fn)
%
% If fn is gzipped, return an ungzipped file and if not, just return the
% input.
% 
% (cc) 2019, sgKIM

[~,~,e1] = myfileparts(fn);
if strcmp(e1(end-1:end),'gz')
  fn_nii = [tempname,e1(1:end-3)];
  system(['gunzip --keep -c > ',fn_nii,' ',fn]);
else
  fn_nii = fn;
end
end