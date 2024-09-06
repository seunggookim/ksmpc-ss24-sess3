function myspm_denan(fname,onlypositive)
nii=load_untouch_nii(fname);
nii.img(isnan(nii.img))=0;
if exist('onlypositive','var')
 nii.img(nii.img<0)=0;
end
save_untouch_nii(nii,fname);
end