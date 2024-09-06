function fname_mat=myspm_coreg_t1w_to_mni(fname_t1w,fname_mni)
% fname_mat=myspm_coreg_t1w_to_mni(fname_t1w,fname_mni)
%
% does skull-stripping and coregistration of a given T1w image to MNI152 space
% * SPM12 must be in the path!
%
% (cc) 2018, sgKIM, solleo@gmail.com

job1=myspm_seg12(struct('fname_t1w',fname_t1w,'mw',zeros(1,6)));
[p1,f1,~]=myfileparts(fname_t1w);
fname_mat=[p1,'/',f1,'_to_MNI.mat'];
if ~exist('fname_mni','var')
fname_fixed=[getenv('FSLDIR'),'/data/standard/MNI152_T1_1mm_brain.nii.gz'];
else
fname_fixed=fname_mni;
end
copyfile(fname_fixed,p1);
gunzip([p1,'/MNI152_T1_1mm_brain.nii.gz']);
M = myspm_coreg_est(job1.fname_out, [p1,'/MNI152_T1_1mm_brain.nii'], fname_mat);
disp('Transform matrix:');
disp(M.T)
disp([' is saved in: ',fname_mat]);

end


function JOB = myspm_seg12(JOB)
% JOB = myspm_seg12(JOB)
%
% JOB requires:
%  .fname_t1w   [string]
%  .isMp2rage  [1x1] default = 0 (if 1 use uncorrected input)
%  .iseastern  [1x1] default = 0 (if 1 use East-Asian template for affine transform)
% (.mw)
% (.naitve)
%
% runs SPM12 segmentation
%
% if isMp2rage=1, then using __original image__ instead of intensity-bias
% corrected one (i.e., m*)
%
% (cc) 2017, sgKIM, solleo@gmail.com

ver=spm('version');
if ~strcmp(ver(1:5),'SPM12')
 error('This function is for SPM12!')
end
if nargin == 0, help(mfilename); return; end

if ~isfield(JOB,'ismp2rage'), JOB.ismp2rage=0; end
if ~isfield(JOB,'iseastern'), JOB.iseastern=0; end
[p2,f2,e2]=myfileparts(JOB.fname_t1w);
dir_tpm=[spm('dir'),'/tpm'];

preproc=[];
preproc.channel.vols = {[JOB.fname_t1w,',1']};
preproc.channel.biasreg = 0.001;
preproc.channel.biasfwhm = 60;
preproc.channel.write = [0 1];
preproc.tissue(1).tpm = {[dir_tpm,'/TPM.nii,1']};
preproc.tissue(1).ngaus = 1;
preproc.tissue(1).native = [1 0];
preproc.tissue(1).warped = [0 1];
preproc.tissue(2).tpm = {[dir_tpm,'/TPM.nii,2']};
preproc.tissue(2).ngaus = 1;
preproc.tissue(2).native = [1 0];
preproc.tissue(2).warped = [0 1];
preproc.tissue(3).tpm = {[dir_tpm,'/TPM.nii,3']};
preproc.tissue(3).ngaus = 2;
preproc.tissue(3).native = [1 0];
preproc.tissue(3).warped = [0 1];
preproc.tissue(4).tpm = {[dir_tpm,'/TPM.nii,4']};
preproc.tissue(4).ngaus = 3;
preproc.tissue(4).native = [0 0];
preproc.tissue(4).warped = [0 0];
preproc.tissue(5).tpm = {[dir_tpm,'/TPM.nii,5']};
preproc.tissue(5).ngaus = 4;
preproc.tissue(5).native = [0 0];
preproc.tissue(5).warped = [0 0];
preproc.tissue(6).tpm = {[dir_tpm,'/TPM.nii,6']};
preproc.tissue(6).ngaus = 2;
preproc.tissue(6).native = [0 0];
preproc.tissue(6).warped = [0 0];
if isfield(JOB,'mw')
 for c=1:6
  preproc.tissue(c).warped(2) = JOB.mw(c);
 end
end
if isfield(JOB,'native')
 for c=1:6
  preproc.tissue(c).native(1) = JOB.native(c);
 end
end
preproc.warp.mrf = 1;
preproc.warp.cleanup = 1;
preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
if JOB.iseastern
 preproc.warp.affreg = 'eastern';
else
 preproc.warp.affreg = 'mni';
end
preproc.warp.fwhm = 0;
preproc.warp.samp = 3;
preproc.warp.write = [0 1];

matlabbatch={};
matlabbatch{1}.spm.spatial.preproc = preproc;

ls(JOB.fname_t1w);
if JOB.ismp2rage
 fname_out = [p2,'/b',f2,e2];
else
 fname_out = [p2,'/bm',f2,e2];
end
if ~exist(fname_out,'file')
 spm_jobman('initcfg');
 spm_jobman('run', matlabbatch);
 Y=cell(1,3);
 for c=1:3
  V = spm_vol_nifti([p2,'/c',num2str(c),f2,e2]);
  [Y{c},XYZ] = spm_read_vols(V);
 end
 if JOB.ismp2rage
  V=spm_vol_nifti([p2,'/',f2,e2]);
 else
  V=spm_vol_nifti([p2,'/m',f2,e2]);
 end
 manat = spm_read_vols(V);
 manat(~(Y{1} | Y{2} | Y{3}))=0;
 V.fname = fname_out;
 V.descrip = 'Skull stripped';
 spm_write_vol(V, manat);
end
disp('> Done:');
ls(fname_out)
JOB.fname_out = fname_out;
end

function out = myspm_coreg_est(fname_moving, fname_fixed, fname_mat)
% T = myspm_coreg_est(fname_moving, fname_fixed, fname_mat)
%
% returns a rigid-body transformation matrix from [fname_moving] to [fname_fixed]
% [fname_mat] can be optionally defined.
%
% (cc) 2017, sgKIM.

if ~nargin,  help myspm_coreg_est;  return; end
[p1,f1,~]=fileparts(fname_moving);
[~,f2,~]=fileparts(fname_fixed);
if ~exist('fname_mat','var')
 fname_mat=[p1,'/coregmat_',f1,'_to_',f2,'.dat'];
end
if ~exist(fname_mat,'file')
 VG=spm_vol(fname_fixed);
 VF=spm_vol(fname_moving);
 x = spm_coreg(VG,VF);
 T = inv(VF.mat\spm_matrix(x(:)')*VG.mat);
 save(fname_mat,'T','-ascii');
else
 T=load(fname_mat,'T','-ascii');
end
out=[];
out.T=T;
out.fname_mat = fname_mat;
end

function [p1,f1,e1]=myfileparts(fname)

[p1,f1,e1] = fileparts(fname);
if isempty(p1) % if not pathname is given
 p1=pwd;
end
if ~strcmp(p1(1),'/') % if the path is relative
 p1=[pwd,'/',p1];
end
if strcmp(e1,'.gz') % if the filename is nii.gz or gii.gz or similar
 e1=[f1(end-3:end),e1];
 f1=f1(1:end-4);
end

end
