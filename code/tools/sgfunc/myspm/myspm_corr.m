function JOB = myspm_corr(JOB)
% JOB = myspm_corr(JOB)
% compute correlation from a BIG file...
%
% JOB reqires:
%  .subjid
%  .dir_base
%  .name_epi
%  .name_voi
%
% (cc) 2015, sgKIM  mailto:solleo@gmail.com


% %% now compute the correlation...
% clear
% load /scr/vatikan4/conmus3/mat/info22.mat
% subjid='AS3T';

subjid   = JOB.subjid;
dir_base = JOB.dir_base;
name_voi = JOB.name_voi;
name_epi = JOB.name_epi;

%dir1=['/scr/vatikan4/conmus3/GLM/',subjid,'/'];
dir1=[dir_base,'/',subjid,'/'];
fname = [dir1,name_epi];
hdr = load_untouch_header_only(fname);
d = hdr.dime.dim(2:5);

%dir2=[dir1,'/corr_STG_right_LMX'];
dir2=[dir1,'/corr_',name_voi];
[~,~]=mkdir(dir2);

%load(['/scr/vatikan4/conmus3/GLM/',subjid,'/glm_FC+FD+BC+BD_RP/VOI_STG_right_LMX_r5mm_1.mat'],'Y');
load([dir_base,subjid,'/glm_FC+FD+BC+BD_RP/VOI_',name_voi,'_r5mm_1.mat'],'Y');
phi=Y;
clear Y
strlength=0;
for z=1:d(3)
nii = load_uns_nii(fname,[], [], [], [], [], z);
nii1 = nii;
nii1.hdr.dime.dim(5) = 1;
nii1.hdr.hist.descrip = ['correlation, z=',pad(z,4),'+',name_voi,'_r5mm'];

if sum(nii.img(:)) == 0
img = nii.img;
else
y = reshape(nii.img, [], d(4))';
r = corr(y, phi);
img = reshape(r', [d(1) d(2)]);
end
nii1.img = img;
fname1=[dir2,'/corr_z',pad(z,4),'.nii'];
str = sprintf('Processing plane: #%d/%d ...\r',z, d(4));
fprintf([repmat('\b',[1,strlength]),str]);
strlength = length(str);
save_untouch_nii(nii1, fname1)
end

unix(['fslmerge -z ',dir2,'/corr_',name_voi,'.nii ',dir2,'/corr_z????.nii']);
unix(['rm -f ',dir2,'/corr_z????.nii'])
