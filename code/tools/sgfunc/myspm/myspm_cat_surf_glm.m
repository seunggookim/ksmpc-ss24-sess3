function JOB = myspm_cat_surf_glm (JOB)
% JOB=myspm_glm (JOB)
%
% This script helps you to set and run GLMs. A result report of
% the 1st-level GLM will be created using myspm_result.m and myspm_graph.m
%
% JOB requires for myspm_glm.m:
% -output directory
%  .dir_base     '1xN' directory for a subdirectory that has SPM.mat
% (.dir_glm)     '1xN' directory to save SPM results
% (.dir_prefix)  '1xN' prefix for GLM directory name
%
% -input files
%  .subjID       'NxM' or {Nx1}
%  .files_query  '1xN' a query to find image filenames with "${subj}",
%                 which will be replaced by .subjID
% or
%  .filenames    {Nsubjx1} (instead of files_query)
%
% (.fwhm_mm)     [1x1|1x3] 3-D smoothing kernel size in mm
% (.masking)     '1xN' filename for an explicit (inclusive) mask
%  .design       '1xN' GLM design: either multiple regression ('mreg')
%                 or one-sample t-test ('t1') or paired t-test ('pt')
%
% -model specification for JOB.design='mreg':
%  .vi.val       [Nsubjx1] a vector of interest
%  .vi.name      'string' a name of interest
%  .vn(c).val    [Nsubjx1] a vector of c-th nuissance variable
%  .vn(c).name   'string' a name of c-th nuissance variable
% or
%  .model        <term> SurfStat term structure that describes a GLM
%  .cidx         [1x1] 1-based index for the contrast of interest
%
% optionally for myspm_result.m:
% (.thres.desc)    '1xN'  'FWE','none', or 'cluster'(default)
% (.thres.alpha)   [1x1]  alpha level (default=0.05)
% (.thres.extent)  [1x1]  extent threshold of clusters in voxels (default=0)
% (.thres.clusterInitAlpha) [1x1] cluster forming height threshold (default=0.001)
% (.thres.clusterInitExtent)  <1x1> cluster forming extent (in voxels) threshold (default=10)
% (.fname_struct)   '1xN' fullpath filename for background anatomical image for orthogonal slices
%                         (defulat='$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz')
% (.titlestr)       {1xNcont} Title text for SPM result report (default={'positive','negative'})
% (.dir_sum)        '1xN' a summary directory into where you want to copy significant results
% (.append)         [1x1] whether to append results into an existing report (default=0)
% (.print)          [1x1] whether to generate results (default=1)
% (.mygraph.x_name) '1xN' x-axis label in a scatterplot (default='x')
% (.mygraph.y_name) '1xN' y-axis label in a scatterplot (default='y')
% (.atlas)          '1xN' atlas to find anatomical names: 'fsl' (default) or 'spm12'
% (.keepResidual)   [1x1] Keeping standardized residuals
%
% Example: one-sample t-test with a covariate
%
% JOB=[];
% JOB.dir_base = '/where/I/want/to/create/dir_glm/';
% JOB.subjID   = {'subj1','subj2'};
% JOB.files_query = '/where/I/have/data/${subj}/fmridata/preproced_epi.nii';
% JOB.design  = 't1';
% JOB.vn.val  = age;
% JOB.vn.name = 'age';
% myspm_glm(JOB)
%
% Results:
%
% (cc) 2015, 2019. sgKIM.  mailto://solleo@gmail.com  https://ggooo.wordpress.com/

if nargin == 0, help(mfilename); return; end
if ~isfield(JOB,'overwrite'), JOB.overwrite=0; end
overwrite=JOB.overwrite;
if ~isfield(JOB,'design'),design='mreg'; else design = JOB.design; end
spm('Defaults','fmri')

matlabbatch={};
matlabbatch{1}.spm.stats.factorial_design.dir = {};
% check data
if isfield(JOB,'fnames')
  fnames = JOB.fnames;
elseif isfield(JOB,'files_query')
  JOB.subjID=fsss_subjID(JOB.subjID);
  fnames=cell(numel(JOB.subjID),1); % this MUST be a column vector <Nx1>
  idx = strfind(JOB.files_query,'${subj}');
  prefix = JOB.files_query(1:idx-1);
  suffix = JOB.files_query(idx+7:end);
  for n=1:numel(JOB.subjID)
    [~,res] = mydir([prefix,JOB.subjID{n},suffix]);
    if isempty(res)
      error(['File not found: ',prefix,JOB.subjID{n},suffix]);
    end
    fnames{n,1}=[res];
  end
elseif isfield(JOB,'filenames')
  for n=1:size(JOB.filenames,1)
    fnames{n,1} = [JOB.filenames{n,1}];
  end
else
  error('You need to specify inputs in JOB.files_query or JOB.filenames');
end
Nsubj=numel(fnames);
for n=1:Nsubj
  ls(fnames{n,1});
end

% check second scans for paired t-test
if strcmpi(design,'pt')
  if isfield(JOB,'files_query2')
    idx = strfind(JOB.files_query2,'${subj}');
    prefix = JOB.files_query2(1:idx-1);
    suffix = JOB.files_query2(idx+7:end);
    for n=1:numel(JOB.subjID)
      [~,res] = mydir([prefix,JOB.subjID{n},suffix]);
      if isempty(res)
        error(['File not found: ',prefix,JOB.subjID{n},suffix]);
      end
      fnames{n,2} = [res];
    end
  elseif size(JOB.filenames,2)
    for n=1:size(JOB.filenames,1)
      fnames{n,2} = [JOB.filenames{n,2}];
      
    end
  else
    error('You need to specify inputs in JOB.files_queryend or JOB.filenames');
  end
  for n=1:Nsubj
    ls(fnames{n,2});
  end
end

%% design specification
switch design
  case 'mreg' % multiple regression
    JOB = myspm_strcNterm(JOB);
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = fnames; % THIS HAS TO BE <Nx1>, not <1xN>!!!
    % variable of interest
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = double(JOB.vi.val);
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = JOB.vi.name;
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 1; % centering
    % including intercept
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
    if isfield(JOB,'nointercept')
      matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 0;
    end
    matlabbatch{1}.spm.stats.factorial_design.multi_cov ...
      = struct('files', {}, 'iCFI', {}, 'iCC', {}); % for spm12
    
  case 't1' % one-sample t-test
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = fnames;
    JOB.vi.name='1';
    
  case 'pt' % paired t-test
    for n=1:size(fnames,1)
      matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(n).scans ...
        = {fnames{n,1};fnames{n,2}};
    end
    matlabbatch{1}.spm.stats.factorial_design.des.pt.gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.pt.ancova = 0;
    JOB.vi.name='1';
end
if isfield(JOB,'vn') % of variable of nuissance
  for c=1:numel(JOB.vn)
    matlabbatch{1}.spm.stats.factorial_design.cov(c).c = double(JOB.vn(c).val);
    matlabbatch{1}.spm.stats.factorial_design.cov(c).cname = JOB.vn(c).name;
    matlabbatch{1}.spm.stats.factorial_design.cov(c).iCFI = 1; % interaction none
    matlabbatch{1}.spm.stats.factorial_design.cov(c).iCC = 1;  % centing overall mean
  end
else
  matlabbatch{1}.spm.stats.factorial_design.cov ...
    = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
end
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
if ~isfield(JOB,'model_desc')
  JOB.model_desc=design;
end
if ~isfield(JOB,'dir_glm')%&&isfield(JOB,'dir_base');
  if ~isfield(JOB,'dir_prefix'), JOB.dir_prefix=''; end
  JOB.dir_glm=fullfile(JOB.dir_base,[JOB.dir_prefix,JOB.model_desc]);
end
[~,~]=mkdir(JOB.dir_glm);
matlabbatch{1}.spm.stats.factorial_design.dir = {JOB.dir_glm};
save([JOB.dir_glm,'/glm_design.mat'], 'matlabbatch');
need2est=1;
if overwrite
  unix(['rm -f ',JOB.dir_glm,'/SPM.mat']);
else
  if exist([JOB.dir_glm,'/SPM.mat'],'file')
    need2est=0;
  end
end
matlabbatch{2}.spm.tools.cat.tools.check_SPM.spmmat{1} = [JOB.dir_glm,'/SPM.mat'];
matlabbatch{2}.spm.tools.cat.tools.check_SPM.check_SPM_cov.do_check_cov.use_unsmoothed_data = 1;
matlabbatch{2}.spm.tools.cat.tools.check_SPM.check_SPM_cov.do_check_cov.adjust_data = 1;
matlabbatch{2}.spm.tools.cat.tools.check_SPM.check_SPM_ortho = 1;
if need2est
  spm_jobman('initcfg')
  spm_jobman('run', matlabbatch)
end
cd(JOB.dir_glm)
spm_print;

%% SURFACE ESTIMATE (set reference mesh correctly (ie, CAT-fsavg-32k/64k) 
swd = JOB.dir_glm;
load([swd,'/SPM.mat'],'SPM')
SPM.swd = swd;
cat_stat_spm(SPM)

%% when required, create a NIFTI file before deleting ResI.*
pwd0=pwd;
cd (JOB.dir_glm)
if isfield(JOB,'keepResidual') && JOB.keepResidual
  % create NIFTI
  setenv('FSLOUTPUTTYPE','NIFTI_GZ')
  system('fslmerge -t ResI ResI_????.nii');
  ls('ResI.nii.gz')
  setenv('FSLOUTPUTTYPE','NIFTI')
end
% delete
system('rm -f ResI_????.*');
cd (pwd0);

%% Contrast estimation
if ~isfield(JOB,'titlestr')
  JOB.titlestr={['+',JOB.vi.name],['-',JOB.vi.name]};
end
if ~isfield(JOB,'cntrstMtx')
  if isfield(JOB,'vn')
    X_cov=zeros(1,numel(JOB.vn));
  else
    X_cov=[];
  end
  if isfield(JOB,'vi')
    if isfield(JOB.vi,'val')
      if isempty( JOB.vi.val )
        X_int=[];
      elseif JOB.vi.val == 1
        X_int=[];
      else
        X_int=0;
      end
    end
  else
    X_int=0;
  end
  switch design
    case 'mreg'
      JOB.cntrstMtx = [X_int X_cov +1; X_int X_cov -1];
    case 't1'
      JOB.cntrstMtx = [+1 X_cov; -1 X_cov];
    case 'pt'
      JOB.cntrstMtx = [
        1 -1 X_cov zeros(1,size(fnames,1));
        -1 1 X_cov zeros(1,size(fnames,1))];
  end
end
JOB.NOREPORT = 1;
myspm_cntrst(JOB);

%% TFCE only the fisrt contrast (for now)

if ~isfield(JOB,'n_perm')
  JOB.n_perm = 5000;
end
if ~isfield(JOB,'is2D')
  JOB.is2D = true;
end


% matlabbatch{1}.spm.tools.tfce_estimate.data = '<UNDEFINED>';
% matlabbatch{1}.spm.tools.tfce_estimate.nproc = 8;
% matlabbatch{1}.spm.tools.tfce_estimate.mask = '';
% matlabbatch{1}.spm.tools.tfce_estimate.conspec.titlestr = '';
% matlabbatch{1}.spm.tools.tfce_estimate.conspec.contrasts = Inf;
% matlabbatch{1}.spm.tools.tfce_estimate.conspec.n_perm = 5000;
% matlabbatch{1}.spm.tools.tfce_estimate.nuisance_method = 2;
% matlabbatch{1}.spm.tools.tfce_estimate.tbss = 0;
% matlabbatch{1}.spm.tools.tfce_estimate.E_weight = 0.5;
% matlabbatch{1}.spm.tools.tfce_estimate.singlethreaded = 0;

matlabbatch={};
for iCnt=1:size(JOB.cntrstMtx,1)
  tfce_estimate = [];
  tfce_estimate.data = {[swd,'/SPM.mat']};
  tfce_estimate.nproc = 8;
  tfce_estimate.mask = '';
  tfce_estimate.conspec.titlestr = JOB.titlestr{iCnt};
  tfce_estimate.conspec.contrasts = iCnt;
  tfce_estimate.conspec.n_perm = JOB.n_perm;
  tfce_estimate.nuisance_method = 2;
  tfce_estimate.tbss = JOB.is2D; % 2D
  tfce_estimate.E_weight = 0.5;
  tfce_estimate.singlethreaded = 0;
  matlabbatch{iCnt}.spm.tools.tfce_estimate = tfce_estimate;
end
spm_jobman('initcfg')
spm_jobman('run', matlabbatch);

%% (convert ps to pdf)
gsbin='/usr/bin/ghostscript';
if ~exist(gsbin,'file')
  gsbin='/opt/local/bin/gs';
end
fname_spm_fig = ['spm_',datestr(now,'yyyymmmdd'),'.ps'];
ps2pdf('psfile',fname_spm_fig, ...
  'pdffile',[fname_spm_fig(1:end-2),'pdf'],...
  'gscommand',gsbin)
delete(fname_spm_fig)

end
