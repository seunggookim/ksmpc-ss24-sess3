function JOB = myspm_ppic (JOB)
% JOB = myspm_ppic (JOB)
%
% JOB requires:
%  .dir_base   = [dir0,mrtID{i}];
%  .dir_voi
% (.dir_psy)
%  .filenames  = {[dir0,mrtID{i},'/swudata.nii']};
%  .voi
%   .name
%   .coord   (mm)
%   .radius  (mm)
%  .ppi
%   .NAMES
%   .CNTRSTMTX  [i j w];
%       first column: i-th condition in SPM.Sess.U(i)
%       second column: j-th parameter? in SPM.Sess.U(i).name{j} ('SOUND', 'SOUNDxLIKE')
%       third column: weight of w
% or
%   .cntrstVec   only contrast weight
%
% for myspm_fmriglm, JOB requires:
%  .TR         (sec)
%  .StimDurSec (sec)
%
%
% shape of the VOI, name, and other options
%
% (cc) 2015, sgKIM.  solleo@gmail.com  https://ggooo.wordpress.com

%% 0. parsing & check inputs

%global overwrite; if isempty(overwrite), overwrite=0; end
if isfield(JOB,'overwrite')&&JOB.overwrite, overwrite=1; else overwrite=0; end
if ~isfield(JOB,'dir_voi')&&isfield(JOB,'dir_psy')
JOB.dir_voi = JOB.dir_psy;
else
JOB.dir_psy = JOB.dir_voi;
end
for n=1:numel(JOB.filenames)
if ~strcmp(JOB.filenames{n}(1),'/') % relative path?
JOB.filenames{n} = fullfile(JOB.dir_base, JOB.filenames{n});
end
fnames{n,1} = [JOB.filenames{n},',1'];
end
try ls(fnames{n,1}(1:end-2));
catch ME
error(['file ',fnames{n,1}(1:end-2),' not found']);
end
JOB.NumSess = numel(fnames);
for j=1:JOB.NumSess
hdr   = load_untouch_header_only(fnames{j}(1:end-2));
NF(j) = hdr.dime.dim(5);
end
JOB.NumFrames=NF;

%% isotropic smoothing
if isfield(JOB,'fwhm')
JOB.fnames = fnames;
JOB = myspm_smooth(JOB);
fnames = JOB.fnames;
end

%% 1. create VOI (physiological factor)

matlabbatch={};
matlabbatch{1}.spm.util.voi.spmmat  = {[JOB.dir_voi,'/SPM.mat']};
matlabbatch{1}.spm.util.voi.adjust  = 0;
matlabbatch{1}.spm.util.voi.session = 1;
matlabbatch{1}.spm.util.voi.name    = [JOB.voi.name,'_r',num2str(JOB.voi.radius),'mm'];
if isfield(JOB.voi,'local') % shift to nearest local maxima in each individual
%   matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat   = {[JOB.dir_voi,'/SPM.mat']};
%   matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = JOB.voi.cntrstVec;
%   matlabbatch{1}.spm.util.voi.roi{1}.spm.conjunction = 1;
%   matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
%   matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.001;
%   matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 0;
%   matlabbatch{1}.spm.util.voi.roi{1}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
%   matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = JOB.voi.coord;
%   matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = JOB.voi.radius;
%   matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.local.spm  = 1;  % this is index of .util.voi.roi
%   matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.local.mask ='';  % and what's this?
%   matlabbatch{1}.spm.util.voi.expression = 'i1&i2';
%   % the updated coordinate is stored in VOI_*_1.mat, xY.xyz
else % fixed coordinate
matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = JOB.voi.coord;
matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = JOB.voi.radius;
matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
matlabbatch{1}.spm.util.voi.expression = 'i1';
end

JOB.fname_voi=[JOB.dir_voi,'/VOI_',JOB.voi.name,'_r',num2str(JOB.voi.radius),'mm_1.mat'];
if ~exist(JOB.fname_voi,'file') || overwrite
spm_jobman('initcfg');
save([JOB.dir_voi,'/mb_voi_pca.mat']);
spm_jobman('run',matlabbatch);
end

%% 2. compute psychological factor (contrast),
%  & product (PPI) of the neural response of physiological factor (deconvoluted BOLD) and psy

num_reg=numel(JOB.ppi.NAMES);
fname_ppi={};
for k=1:num_reg
% first column:  i in SPM.Sess.U(i)
% second column: j in SPM.Sess.U(i).name{j}
% third column: weight
if ~isfield(JOB,'dir_ppi')
dir_ppi = [JOB.dir_base,'/ppi_',JOB.ppi.dirname];
end
[~,~] = mkdir(dir_ppi);
matlabbatch={};
matlabbatch{1}.spm.stats.ppi.spmmat = {[JOB.dir_psy,'/SPM.mat']};
matlabbatch{1}.spm.stats.ppi.type.ppi.voi = {JOB.fname_voi};
JOB.ppi.CNTRSTMTX{k}
matlabbatch{1}.spm.stats.ppi.type.ppi.u   = JOB.ppi.CNTRSTMTX{k};
matlabbatch{1}.spm.stats.ppi.name = JOB.ppi.NAMES{k}; %'STG-LxBCFD';
matlabbatch{1}.spm.stats.ppi.disp = 0; % not very informative
fname_ppi{k}=[dir_ppi,'/PPI_',JOB.ppi.NAMES{k},'.mat'];
if ~exist(fname_ppi{k},'file') || overwrite
spm_jobman('initcfg')
save([dir_ppi,'/mb_ppi_deconv.mat'],'matlabbatch');
spm_jobman('run',matlabbatch);
fname_output=[JOB.dir_psy,'/PPI_',JOB.ppi.NAMES{k},'.mat'];
copyfile(fname_output,fname_ppi{k})
end
end
JOB.fname_ppi = fname_ppi;

%% 3. combine regressors [phy, psy1*phy, psy2*phy, psy3*phy, ... ]

load(JOB.fname_ppi{1}, 'PPI');
JOB.reg(1).name = 'PHY';
JOB.reg(1).val  = PPI.Y;
JOB.cntrstMtx = 0;
num_reg  = numel(JOB.fname_ppi);
PSYNAMES = JOB.ppi.NAMES;
i=2;
for k=1:num_reg
load(JOB.fname_ppi{k}, 'PPI');
JOB.reg(i).name = PSYNAMES{k};
JOB.reg(i).val  = PPI.P;
JOB.reg(i+1).name = [PSYNAMES{k},'xPHY'];
JOB.reg(i+1).val  = PPI.ppi;
i=i+2;
JOB.cntrstMtx=[JOB.cntrstMtx 0 JOB.PPICntrstVec(k)];
end
JOB.cntrstMtx=[JOB.cntrstMtx 0]; % for interceptor
JOB.cntrstMtx=[JOB.cntrstMtx; -JOB.cntrstMtx];

%% 4. run fmriglm!

JOB.dir_glm = dir_ppi;
JOB = myspm_fmriglm(JOB);
