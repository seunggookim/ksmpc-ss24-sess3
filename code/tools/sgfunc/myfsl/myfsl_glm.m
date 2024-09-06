function JOB = myfsl_glm (JOB)
% JOB = myfsl_glm (JOB)
%
% sets glm and runs randomization for VBM and TBSS (fname_input with
% "skeletonised")
%
% JOB requires:
%  .model        [term] SurfStat 'term' format of a GLM
%  .cidx         [1x1] the index of a regressor of interest in the model
% (.numPerm)     [1x1] the number of randomization
%  .dir_base     'Nx1' a directory that encloses a model's subdirectory
%  .fname_input  'Nx1' fullpath filename of "all_*"
%  .fname_mask   'Nx1' fullpath filename of a mask
% (.glm_prefix)
%
% (cc) 2015, sgKIM.  solleo@gmail.com   https://ggooo.wordpress.com
if ~isfield(JOB,'glm_prefix'), JOB.glm_prefix='';  end

if ~isfield(JOB, 'model_desc') && ~isfield(JOB,'COI')
 [JOB.model_desc,JOB.COI] = fsss_model_desc(JOB.model, JOB.cidx);
end

% create output directory
if ~isfield(JOB,'dir_glm')
 JOB.dir_glm = [JOB.dir_base,'/GLM/',JOB.meas,'_',JOB.glm_prefix,JOB.model_desc,'/'];
end
[~,~] = mkdir(JOB.dir_glm);

% create design and contrasts files
JOB = myfsl_designmat(JOB);
JOB.output_prefix = fullfile(JOB.dir_glm,JOB.vi);

% 10K randomization as a default
if ~isfield(JOB,'numPerm'),  JOB.numPerm=10*1000; end

% is it TBSS?
if strfind(JOB.fname_input,'skeletonised')
 option_2d = ' --T2 '; % for tbss
else
 option_2d = ' -T ';
end

% run glm until it's done?
rerun=[1 1];
while ~~sum(rerun)
 for s=1:2
  fname = [JOB.dir_glm,'/',JOB.COI,'_tfce_corrp_tstat',num2str(s),'.nii.gz'];
  if exist(fname,'file')
   rerun(s) = 0;
  end
 end
 %[~,qid] = myunix(['cd ',JOB.dir_glm,'; randomise_parallel -i ',JOB.fname_input, ...
 [~,qid] = myunix(['cd ',JOB.dir_glm,'; randomise -i ',JOB.fname_input, ...
  ' -o ',JOB.output_prefix,' -m ',JOB.fname_mask, ...
  ' -d ',JOB.fname_mat,' -t ',JOB.fname_con,' -n ',num2str(JOB.numPerm),...
  option_2d],1);
 if strfind(qid,'ERROR')
  qid
  error('Error occured in submitting to condor!');
 end
 qid = qid(end-4:end-1);
 disp(['# GLM:',JOB.meas,':',JOB.model_desc,':',JOB.COI,' submitted qid=',qid]);
 JOB.qid = qid;
 if isfield(JOB,'getItDone')
  if ~~sum(rerun)
   myunix(['waitForCONDORJobs.sh 60 ',JOB.qid]);
  end
 else
  rerun=0;
 end
end

end
