function nii = load_untouch_nii_like(fname, ref)
% nii = load_untouch_nii_like(fname, ref)
% (cc) sgKIM

fname_tmp1 = [tempname,'.nii'];
[a,b]=unix(['FSLOUTPUTTYPE=NIFTI; fslroi ',fname,' ',fname_tmp1,' 0 1']);
fname_tmp = [tempname,'.nii'];
[a,b]=unix(['mri_convert --resample_type nearest --like ',ref,' ',fname_tmp1,' ',fname_tmp]);
if a, disp(b); end
nii = load_uns_nii(fname_tmp);
[~,~]=unix(['rm -f ',fname_tmp1]);
[~,~]=unix(['rm -f ',fname_tmp]);
end
