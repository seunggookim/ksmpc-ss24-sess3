function JOB = myspm_deconv (JOB)
%{
JOB = myspm_deconv (JOB)

requires:
JOB
.dir_base
.dir_phi  (default: dir_base)

.phi
.name
.fname_voi | .fname_mask | .cluster | .sphere (either one of them)
.fname_voi (for already computed VOI .mat file)
or
.fname_mask (for a mask image .nii or .img file)
or
.cluster (for a cluster ID maps (sigclus*.nii) created from myspm_fmriglm)
.index
.sign (either '+' or '-')
or
.sphere (spherical VOI)
.centre (mm)
.radius (mm)

.psy
.name
.cntrstVec (psychological contrast as in GLM)

(.ppi.name) default: [phi.name,'x',psy.name]

OUTPUT:

Time of the first sample is SPM.xBF.T0-th mirco-time sample
If dt is 0.0625 (or 1/16) sec as usual, the T(0) is 0.5 sec.

See:
RT=SPM.xY.RT
dt=SPM.xBF.dt
fMRI_T0=SPM.xBF.T0

(cc) 2015, 2016, sgKIM.  solleo@gmail.com  https://ggooo.wordpress.com
%}

%% 0. parsing & check inputs

global overwrite; if isempty(overwrite), overwrite=0; end
spm('Defaults','fmri');
spm_jobman('initcfg');
for n=1:numel(JOB.filenames)
if ~strcmp(JOB.filenames{n}(1),'/') % relative path?
JOB.filenames{n} = fullfile(JOB.dir_base, JOB.filenames{n});
end
fnames{n,1} = [JOB.filenames{n},',1'];
end
try ls(fnames{n,1}(1:end-2));
catch ME
fnames
error(['file ',fnames{n,1}(1:end-2),' not found']);
end
JOB.NumSess = numel(fnames);
for j=1:JOB.NumSess
hdr   = load_untouch_header_only(fnames{j}(1:end-2));
NF(j) = hdr.dime.dim(5);
end
JOB.NumFrames=NF;

%if ~isfield(JOB,'dir_phi'), JOB.dir_phi = JOB.dir_base; end
% if ~isfield(JOB,'dir_psy'), JOB.dir_psy = JOB.dir_base; end
[p1,~,~]=fileparts(JOB.dir_glm);
JOB.dir_base=p1;
JOB.dir_phi=JOB.dir_glm;

%% isotropic smoothing
if isfield(JOB,'fwhm')
JOB.fnames = fnames;
JOB = myspm_smooth(JOB);
fnames = JOB.fnames;
end

%% 1. create VOI (for physiological vector)
spm('Defaults','fmri');
spm_jobman('initcfg');

matlabbatch={};
matlabbatch{1}.spm.util.voi.spmmat  = {[JOB.dir_phi,'/SPM.mat']};
matlabbatch{1}.spm.util.voi.adjust  = 0;
matlabbatch{1}.spm.util.voi.session = 1;
if ~isfield(JOB,'fname_voi')
if isfield(JOB.phi,'fname_mask') % just a file
matlabbatch{1}.spm.util.voi.name = JOB.phi.name;
matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {[JOB.phi.fname_mask,',1']};
matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0;
matlabbatch{1}.spm.util.voi.expression = 'i1';
JOB.fname_voi=[JOB.dir_phi,'/VOI_',JOB.phi.name,'_1.mat'];
elseif isfield(JOB.phi,'cluster') % create a mask file for i-th contrast, j-th cluster
matlabbatch{1}.spm.util.voi.name = ['clus_',JOB.phi.name];
if ~isfield(JOB.phi.cluster,'dir'); JOB.phi.cluster.dir=JOB.dir_phi; end
fname_sig = [JOB.phi.cluster.dir,'/sigclus_',JOB.phi.cluster.sign,'1.nii'];
nii = load_uns_nii(fname_sig);
cluster_indices=sort(unique(nii.img(:)));
cluster_indices(1)=[];
nii.img = uint8(nii.img == cluster_indices(JOB.phi.cluster.index));
nii.hdr.dime.datatype=2;
fname_cluster = [JOB.dir_phi,'/clusmask_',JOB.phi.name,'.nii'];
save_untouch_nii(nii, fname_cluster);
clear nii
matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {[fname_cluster,',1']};
matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0;
matlabbatch{1}.spm.util.voi.expression = 'i1';
JOB.fname_voi = [JOB.dir_phi,'/', ...
'VOI_',matlabbatch{1}.spm.util.voi.name,'_1.mat'];
elseif isfield(JOB.phi,'sphere') % sphere
matlabbatch{1}.spm.util.voi.name = [JOB.phi.name,'_r',num2str(JOB.phi.sphere.radius),'mm'];
matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = JOB.phi.sphere.centre;
matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = JOB.phi.sphere.radius;
matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
matlabbatch{1}.spm.util.voi.expression = 'i1';
JOB.fname_voi=[JOB.dir_phi,'/', ...
'VOI_',JOB.phi.name,'_r',num2str(JOB.phi.sphere.radius),'mm_1.mat'];
end
% % in case of Y is not residual
%  if isfield(JOB,'voiadjust') && JOB.phiadjust
%   matlabbatch{1}.spm.util.voi.adjust = NaN; % This seems like taking eig1 from residuals... is it a right way to do?
%  end
end
if ~exist(JOB.fname_voi,'file') || overwrite
save([JOB.dir_phi,'/myspm_ppi_1_create_voi.mat'],'matlabbatch');
spm_jobman('run',matlabbatch);
end

%% 2. do simple deconvolution (sd)
matlabbatch={};
matlabbatch{1}.spm.stats.ppi.spmmat = {[JOB.dir_phi,'/SPM.mat']};
matlabbatch{1}.spm.stats.ppi.type.sd.voi = {JOB.fname_voi};
matlabbatch{1}.spm.stats.ppi.name=JOB.phi.name;
matlabbatch{1}.spm.stats.ppi.disp = 0; % not very informative
JOB.fname_ppi=[JOB.dir_phi,'/PPI_',JOB.phi.name,'.mat'];

spm_jobman('initcfg')
if ~exist(JOB.fname_ppi,'file') || overwrite
save([JOB.dir_phi,'/myspm_ppi_2_construct_ppi.mat'],'matlabbatch');
spm_jobman('run',matlabbatch);
%copyfile([JOB.dir_phi,'/PPI_',JOB.ppi.name,'.mat'],dir_ppi);
end

end
