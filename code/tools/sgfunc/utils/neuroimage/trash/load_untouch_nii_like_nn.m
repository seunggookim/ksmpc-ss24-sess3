function nii = load_untouch_nii_like_nn(fname,ref)
% %fname_tmp = ['/tmp/',num2str(randi(100)),'.nii'];
% [p1,n1,e1]=fileparts(fname);
% fname_intp=[p1,'/',n1,'_1mm',e1];
% if ~exist(fname_intp,'file')
% myunix(['mri_convert_nn ',ref,' ',fname,' ',fname_intp],1);
% end
% nii = load_uns_nii(fname_intp);
[p1,f1,e1]=myfileparts(fname);
fn1=[tempname,e1];
copyfile(fname,fn1);
if strcmp(e1,'.nii.gz')
  gunzip(fn1);
  fn1=fn1(1:end-3);
end

[p2,f2,e2]=myfileparts(ref);
fn2=[tempname,e2];
copyfile(ref,fn2)
if strcmp(e2,'.nii.gz')
  gunzip(fn2);
  fn2=fn2(1:end-3);
end

job=myspm_reslice(struct('fname_src',fn1, 'fname_ref',fn2, 'interp',0));
nii= load_untouch_nii(job.fname_out);


end
