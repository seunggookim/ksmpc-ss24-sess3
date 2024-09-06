function JOB = myspm_rsfc_seed(JOB)
% JOB = myspm_rsfc_seed(JOB)
%
% Simple and crisp correlation computing from 4-D fMRI data!
%
% JOB requires:
%  .fname_data
%  .fname_roi
%  .fname_out

% data
fname_data = JOB.fname_data;
ls(fname_data);
hdr=load_untouch_header_only(JOB.fname_data);
d=hdr.dime.dim(2:5);
% roi or u1
if isfield(JOB,'fname_roi') && ~isfield(JOB,'fname_u1')
 fname_roi = JOB.fname_roi;
 ls(fname_roi);
 roi = load_uns_nii(JOB.fname_roi);
 Y = spm_get_data(JOB.fname_data,find3(~~roi.img)');
 u1=pc1(Y);
 [p1,n1,~]=fileparts(JOB.fname_roi);
 JOB.fname_u1=[p1,'/',n1,'.mat'];
 disp(['[',subjid,']: saving PC1: ',JOB.fname_u1]);
 save(JOB.fname_u1,'u1','-ascii');
elseif ~isfield(JOB,'fname_roi') && isfield(JOB,'fname_u1')
 ls(JOB.fname_u1);
 u1=load(JOB.fname_u1);
else
 error('.fname_roi or .fname_u1 field is required!');
end
% output filename
fname_out = JOB.fname_out;
[p1,~,~] = fileparts(fname_out);
[~,~]=mkdir(p1);
% if already done, just exit
if exist(JOB.fname_out,'file') && ~(isfield(JOB,'overwrite') && JOB.overwrite)
 return
end
% computing each slice
R=zeros(d(1:3));
disp(['[',subjid,']: ',num2str(d(3)),' slices to compute']);
for z=1:d(3)
 nii = load_uns_nii(JOB.fname_data, 1:d(4), [], [], [], [], z);
 y = reshape(squeeze(nii.img),[],d(4))';
 r = corr(y, reshape(u1,d(4),[]));
 corrslice = reshape(r,d(1:2));
 R(:,:,z)  = corrslice;
end
% save as .nii
nii0=nii;
nii.hdr.dime.dim=[3 d(1:3) 1 1 1 1];
nii.hdr.dime.datatype=16;
[p1,f1,~]=fileparts(JOB.fname_out);
nii.fileprefix=[p1,'/',f1];
nii.img = R;
%nii.hdr.hist.descrip=['conmus3: ',nii.hdr.hist.descrip,': correlated with inferior colliculus'];
save_untouch_nii(nii, JOB.fname_out)
disp(['[',JOB.subjid,']: ']);

end
