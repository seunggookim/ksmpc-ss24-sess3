function job = myspm_glm(job)
% job = myspm_glm(job)
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
%  .subjID       [NxM] or {Nx1}
%  .files_query  'Nx1' a query to find image filenames with "${subj}",
%                 which will be replaced by .subjID
% or
%  .filenames    {Nsubjx1} (instead of files_query) for 'mreg'
% or
%  .filenames    {NsubjxNcond} for 'flex'
% 
%
% (.masking)     '1xN' filename for an explicit (inclusive) mask
%  .design       '1xN' GLM design: either multiple regression ('mreg')
%                 or one-sample t-test ('t1') or paired t-test ('pt')
%                 or flexible factorial
%
% -factor specification for 'flex'
%  .factors      (1xNcond)  .name (.dept) (.variance) (.gmsca) (.ancova)
%  .conds        {1xNsub} of which element is [NimagesPerSubjectsxNcond]
%  .effects      {1xNeffects} factor index, [1] for main effect of factor 1
%                [1 2] for interaction between factor 1 and factor2
%
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
% -for myspm_cntrst.m:
% (.cntrstMtx)
% (.titlestr)
% (.effectOfInterest)
% (.FcntrstMtx)
% (.Ftitlestr)
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
% (cc) 2015. sgKIM.  mailto://solleo@gmail.com  https://ggooo.wordpress.com/

if nargin == 0, help(mfilename); return; end
if ~isfield(job,'overwrite'), job.overwrite=0; end
overwrite=job.overwrite;
if ~isfield(job,'design'),design = 'mreg'; else, design = job.design; end
spm('Defaults','fmri')

matlabbatch={};
matlabbatch{1}.spm.stats.factorial_design.dir = {};
% check data
if isfield(job,'fnames')
  fnames = job.fnames;
elseif isfield(job,'files_query')
  job.subjID=fsss_subjID(job.subjID);
  fnames=cell(numel(job.subjID),1); % this MUST be a column vector <Nx1>
  idx = strfind(job.files_query,'${subj}');
  prefix = job.files_query(1:idx-1);
  suffix = job.files_query(idx+7:end);
  for n=1:numel(job.subjID)
    [~,res] = mydir([prefix,job.subjID{n},suffix]);
    if isempty(res)
      error(['File not found: ',prefix,job.subjID{n},suffix]);
    end
    fnames{n,1}=[res,',1'];
  end
elseif isfield(job,'filenames')
  for n=1:size(job.filenames,1)
    fnames{n,1} = [job.filenames{n,1},',1'];
  end
else
  error('You need to specify inputs in JOB.files_query or JOB.filenames');
end
Nsubj=numel(fnames);
for n=1:Nsubj
  ls(fnames{n,1}(1:end-2));
end

% check second scans for paired t-test
if strcmpi(design,'pt')
  if isfield(job,'files_query2')
    idx = strfind(job.files_query2,'${subj}');
    prefix = job.files_query2(1:idx-1);
    suffix = job.files_query2(idx+7:end);
    for n=1:numel(job.subjID)
      [~,res] = mydir([prefix,job.subjID{n},suffix]);
      if isempty(res)
        error(['File not found: ',prefix,job.subjID{n},suffix]);
      end
      fnames{n,2} = [res,',1'];
    end
  elseif size(job.filenames,2)
    for n=1:size(job.filenames,1)
      fnames{n,2} = [job.filenames{n,2},',1'];
    end
  else
    error('You need to specify inputs in JOB.files_queryend or JOB.filenames');
  end
  for n=1:Nsubj
    ls(fnames{n,2}(1:end-2));
  end
end

%% design specification
switch design
  case 'mreg' % multiple regression
    job = myspm_strcNterm(job);
    des = [];
    des.mreg.scans = fnames; % THIS HAS TO BE <Nx1>, not <1xN>!!!
    % variable of interest
    des.mreg.mcov.c = double(job.vi.val);
    des.mreg.mcov.cname = job.vi.name;
    des.mreg.mcov.iCC = 1; % centering
    % including intercept
    des.mreg.incint = 1;
    if isfield(job,'nointercept')
      des.mreg.incint = 0;
    end
    matlabbatch{1}.spm.stats.factorial_design.des = des;    
    matlabbatch{1}.spm.stats.factorial_design.multi_cov ...
      = struct('files', {}, 'iCFI', {}, 'iCC', {}); % for spm12
    
  case 't1' % one-sample t-test
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = fnames;
    job.vi.name = '1';
    job.cntrstMtx = [1; -1];
    
  case 'pt' % paired t-test
    for n=1:size(fnames,1)
      matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(n).scans ...
        = {fnames{n,1};fnames{n,2}};
    end
    matlabbatch{1}.spm.stats.factorial_design.des.pt.gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.pt.ancova = 0;
    job.vi.name='1';
    
  case 'flex' % flexible factorial
    des = [];
    des.fblock.fac = job.factors;
    des.fblock.fsuball.fsubject = [];
    for isub = 1:size(job.filenames,1)
      des.fblock.fsuball.fsubject(isub).scans = reshape(...
        job.filenames(isub,1:size(job.conds{isub},1)),[],1);
      des.fblock.fsuball.fsubject(isub).conds = job.conds{isub};
    end
    des.fblock.maininters = {};
    for ieff = 1:numel(job.effects)
      if numel(job.effects{ieff}) == 1
        des.fblock.maininters{ieff}.fmain.fnum = job.effects{ieff};
      else
        des.fblock.maininters{ieff}.inter.fnums = job.effects{ieff}(:);
      end
    end
    matlabbatch{1}.spm.stats.factorial_design.des = des;
end
if isfield(job,'vn') % of variable of nuissance
  for c=1:numel(job.vn)
    matlabbatch{1}.spm.stats.factorial_design.cov(c).c = double(job.vn(c).val);
    matlabbatch{1}.spm.stats.factorial_design.cov(c).cname = job.vn(c).name;
    matlabbatch{1}.spm.stats.factorial_design.cov(c).iCFI = 1; % interaction none
    matlabbatch{1}.spm.stats.factorial_design.cov(c).iCC = 1;  % centing overall mean
  end
else
  matlabbatch{1}.spm.stats.factorial_design.cov ...
    = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
end
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
if ~isfield(job,'masking')
  matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
else
  matlabbatch{1}.spm.stats.factorial_design.masking.em = {[job.masking,',1']};
end
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
if ~isfield(job,'model_desc')
  job.model_desc = design;
end
if ~isfield(job,'dir_glm')%&&isfield(JOB,'dir_base');
  if ~isfield(job,'dir_prefix'), job.dir_prefix = ''; end
  job.dir_glm = fullfile(job.dir_base,[job.dir_prefix,job.model_desc]);
end
[~,~] = mkdir(job.dir_glm);
matlabbatch{1}.spm.stats.factorial_design.dir = {job.dir_glm};
save([job.dir_glm,'/glm_design.mat'], 'matlabbatch');
need2est = 1;
if overwrite
  unix(['rm -f ',job.dir_glm,'/SPM.mat']);
else
  if exist([job.dir_glm,'/SPM.mat'],'file')
    need2est = 0;
  end
end
if need2est
  spm_jobman('initcfg')
  spm_jobman('run', matlabbatch)
end
cd(job.dir_glm)
spm_print;

%% Coefficeint estimation
matlabbatch = {};
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[job.dir_glm,'/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
save([job.dir_glm,'/glm_estimation.mat'], 'matlabbatch');
if overwrite || ~exist([job.dir_glm,'/beta_0001.nii'],'file')
  spm_jobman('run', matlabbatch)
end


%% when required, create a NIFTI file before deleting ResI.*
pwd0 = pwd;
cd (job.dir_glm)
if isfield(job,'keepResidual') && job.keepResidual
  % create NIFTI
  setenv('FSLOUTPUTTYPE','NIFTI_GZ')
  myunix('fslmerge -t ResI ResI_????.nii');
  ls('ResI.nii.gz')
  setenv('FSLOUTPUTTYPE','NIFTI')
end
% delete
unix('rm -f ResI_????.nii');
cd (pwd0);


%% 4. Design review & Orthogonality check
% set figure filename
if ~isfield(job,'fname_spm_fig')
  today = datestr(now,'yyyymmmdd');
  job.fname_spm_fig = fullfile(job.dir_glm,[job.model_desc,'_spm_',today,'.ps']);
  job.fname_spm_fig = strrep(job.fname_spm_fig,'>','-gt-');
  job.fname_spm_fig = strrep(job.fname_spm_fig,'<','-lt-');
end
% delete previous figure files today (rewriting)
if exist(job.fname_spm_fig,'file')
  delete(job.fname_spm_fig);
end

% review the design matrix and save it:
toDisplay = {'matrix','orth','covariance'};
spm_figure('clear')
for j=1:numel(toDisplay)
  matlabbatch = {};
  matlabbatch{1}.spm.stats.review.spmmat = {[job.dir_glm,'/SPM.mat']};
  matlabbatch{1}.spm.stats.review.display.(toDisplay{j})= 1;
  matlabbatch{1}.spm.stats.review.print = false;
  spm_jobman('run', matlabbatch);
  spm_print(job.fname_spm_fig);
end


%% Now create result reports
if ~isfield(job,'thres')
  job.thres.desc  = 'cluster';
  job.thres.alpha = 0.05;
end
if ~(isfield(job,'NOCNTRST') && job.NOCNTRST)
  myspm_cntrst (job);
end

end
