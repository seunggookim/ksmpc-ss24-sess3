function myfsl_smooth_within_mask(JOB)
% myfsl_smooth_within_mask(JOB)
%
% JOB requires:
%  .fname_input
%  .fwhm

sigma= num2str((JOB.fwhm/(2*sqrt(2*log(2)))));
fname_input=JOB.fname_input;
[p1,n1,e1]=fileparts_gz(fname_input);
if isempty(p1), p1=pwd; end
fname_output = [p1,'/aniso_s',num2str(JOB.fwhm),'_',n1,e1];

cmd0='FSLOUTPUTTYPE=NIFTI; fslmaths ';
unix([cmd0,fname_input,' -bin -s ',sigma,' /tmp/binsmooth.nii']);
bs=load_uns_nii('/tmp/binsmooth.nii');
unix([cmd0,fname_input,' -s ',sigma,' /tmp/smooth.nii']);
sm=load_uns_nii('/tmp/smooth.nii');

img=bs.img./sm.img;
img(isnan(img))=0;
img(isinf(img))=0;
M=max(sm.img(:));
m=min(sm.img(:));
img(img>M)=M;
img(img<m)=m;

nii=sm;
nii.img = img;
save_untouch_nii(nii, fname_output);

end
