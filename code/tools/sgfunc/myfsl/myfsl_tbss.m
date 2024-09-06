function JOB = myfsl_tbss (JOB)
% JOB = fsss_glm(JOB)
%
% JOB requires:
%  .subjID       [Nsubjx1]|{Nsubjx1}
%  .model        <term> SurfStat model term
%  .cdix         [1x1] index of term to test
%  .dir_input    (where you have ${dir_input}/${SUBJID}/*${input_prefix}*FA*.nii.gz files
%  .input_prefix 'Nx1'
% or
%  .fname_input {NmeasxNsubj}
%  .dir_base    (default:pwd)
%  .glm_prefix
% (.alpha)       [1x1] alpha level (default: 0.05)
%
% (cc) 2015, sgKIM.   solleo@gmail.com

isPar=1;
if ~nargin,   help myfsl_tbss;  return;  end
if ~isfield(JOB,'alpha'),   JOB.alpha=0.05;  end
subjID= fsss_subjID(JOB.subjID);
N = numel(subjID);
if ~isfield(JOB,'glm_prefix'),  JOB.glm_prefix=''; end
setenv('FSLOUTPUTTYPE','NIFTI_GZ');
if isfield(JOB,'dir_base'),  dir_base=JOB.dir_base;
else  dir_base=pwd;  end
dir_log = [dir_base,'/job_logs'];
[~,~]=mkdir(dir_log);
cd(dir_log);
dir_figs  = [dir_base ,'/fig_overview'];
[~,~]=mkdir(dir_figs);

%% A. First coregst all inputs
if isfield(JOB,'measures')
 measures = {'FA',JOB.measures{:}};
else
 measures = {'FA','L1','L2','L3','MD'};
end
if ~exist([dir_base,'/inputs/',subjID{N},'/r',measures{end},'.nii'],'file')
 for n=1:N
  subjid = subjID{n};
  fixed = '/home/raid2/skim/FMRIB58_FA_1mm.nii';
  src = [dir_base,'/inputs/',subjid,'/FA.nii'];
  job1=struct('name_fixed',fixed, 'name_moving',src, 'interp',0, 'prefix','r');
  for k=2:numel(measures)
   job1.name_others{k-1,1} = [dir_base,'/inputs/',subjid,'/',measures{k},'.nii'];
  end
  myspm_coreg(job1);
  unix(['cd ',dir_base,'/inputs/',subjid,'; FSLOUTPUTTYPE=NIFTI; ', ...
   'fslmaths rL2 -add rL3 -div 2 rL23']);
 end
end
%% B. do all subject with FA first
if ~exist([dir_base,'/stats/mean_FA_skeleton_mask.nii.gz'],'file')
 %% tbss step #1: seting directory/files (EBG; Excluding BackGround from mask)
 if ~exist([dir_base,'/origdata/FA',subjID{N},'.nii'],'file')
  for n=1:N
   subjid = subjID{n};
   src = [dir_base,'/inputs/',subjid,'/rFA.nii'];
   trg = [dir_base,'/',subjid,'.nii'];
   unix(['ln -sf ',src,' ',trg]);
  end
  unix(['cd ',dir_base,'; tbss_1_preproc_EBG *.nii']);
  unix(['cd ',dir_base,'/FA/; slicesdir *_FA.nii.gz; ', ...
   ' mv ',dir_base,'/FA/slicesdir ',dir_figs,'/FA_native']);
 end
 %% tbss step #2: FLIRT and FNIRT
 if ~exist([dir_base,'/FA/',subjID{N},'_FA_to_target_warp.nii.gz'],'file')
  if isPar
   [~,qid]=unix(['cd ',dir_base,'; ',...
    ' tbss_2_reg -T -flirt " -searchrx -15 15 -searchry -15 15 -searchrz -15 15 "']);
   qid = qid(end-4:end-1);
   disp(['# Running FNIRT... qid=',qid]);
   myunix(['waitForCONDORJobs.sh 300 ',qid]);
  else
   unix(['cd ',dir_base,'; tbss_2_reg_aff_serial -T -flirt " -coarsesearch 30 "']);
  end
 end
 %% tbss step #3: applying warping
 if ~exist([dir_figs,'/FA_fnirted'],'dir');
  unix(['cd ',dir_base,'; tbss_3_postreg -S;']);
  % for visual inspection
  unix(['cd ',dir_base,'/FA/; slicesdir *_FA_to_target.nii.gz']);
  myunix(['mv ',dir_base,'/FA/slicesdir ',dir_figs,'/FA_fnirted'],1);
 end
 %% tbss step #4: skeletonizing
 if ~exist([dir_base,'/FA/mean_FA_skeleton_mask.nii.gz'],'file')
  unix(['cd ',dir_base,';. tbss_4_prestats 0.2;']);
  % print skeleton mask overlaid on mean FA template
  myunix(['slices ',dir_base,'/stats/mean_FA.nii.gz ',...
   dir_base,'/stats/mean_FA_skeleton_mask.nii.gz ',...
   ' -o ',dir_figs,'/mean_FA_skeleton_mask_OVER_mean_FA.gif'],1);
 end
end

%% B. non-FA indices
othermeas = measures(2:end);
for k=1:numel(othermeas)
 meas = othermeas{k};
 if ~exist([dir_base,'/stats/all_',meas,'_skeletonised.nii.gz'],'file')
  [~,~]=mkdir([dir_base,'/',meas]);
  for n=1:N
   subjid = subjID{n};
   src = [dir_base,'/inputs/',subjid,'/r',meas,'.nii'];
   trg = [dir_base,'/',meas,'/',subjid,'.nii'];
   unix(['ln -sf ',src,' ',trg]);
  end
  unix(['cd ',dir_base,'/',meas,';slicesdir *;', ...
   ' mv slicesdir ',dir_figs,'/',meas,'_native']);
  
  unix(['cd ',dir_base,'; tbss_non_FA ',meas]);
  unix(['cd ',dir_base,'/FA; rm -rf slicesdir;', ...
   ' slicesdir ????_to_target_',meas,'.nii*;', ...
   ' mv slicesdir ',dir_figs,'/',meas,'_fnirted']);
 end
end


%% C. GLM
numMeas = numel(measures);
rerun=ones(1,numMeas);
while ~~sum(rerun)
 for k=1:numMeas
  meas = measures{k};
  [JOB.model_desc,JOB.COI] = fsss_model_desc(JOB.model,JOB.cidx);
  % get links from "stats" directory to "GLM" directories
  dir_glm = [dir_base,'/GLM/',meas,'_',JOB.glm_prefix,JOB.model_desc,'/'];
  for s=1:2
   X = char(JOB.model);
   fname = [dir_glm,JOB.COI,'_tfce_corrp_tstat',num2str(s),'.nii.gz'];
   if exist(fname,'file')
    rerun(k) = 0;
   end
  end
 end
 for k=1:numMeas
  meas = measures{k};
  JOB.model_desc = fsss_model_desc(JOB.model,JOB.cidx);
  % get links from "stats" directory to "GLM" directories
  dir_glm = [dir_base,'/GLM/',meas,'_',JOB.glm_prefix,JOB.model_desc,'/'];
  [~,~] = mkdir(dir_glm);
  cd(dir_glm);
  fnames=dir([dir_base,'/stats/*',meas,'*']); % all_FA.nii.gz, mean_FA.nii.gz, ...
  for j=1:numel(fnames)
   if ~fnames(j).isdir
    unix(['ln -sf ',dir_base,'/stats/',fnames(j).name,' ',dir_glm,'.']);
   end
  end
  JOB.fname_input = [dir_glm,'/all_',meas,'_skeletonised.nii.gz'];
  ls(JOB.fname_input);
  
  JOB.fname_mask  = [dir_glm,'/mean_FA_skeleton_mask.nii.gz'];
  unix(['ln -sf ',dir_base,'/stats/mean_FA_skeleton_mask.nii.gz ',JOB.fname_mask]);
  ls(JOB.fname_mask);
  
  JOB.dir_glm = dir_glm;
  if rerun(k)
   JOB = myfsl_glm(JOB);
  end
 end
 if ~~sum(rerun)
  myunix(['waitForCONDORJobs.sh 60 ',JOB.qid]);
 end
end
%% need to confirm everything was done


end


