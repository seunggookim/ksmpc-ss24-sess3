function myspm_par(funcname, JOBs, num_workers)
% myspm_par(funcname, JOBs, num_workers)
% (cc) 2019, sgKIM.

func = str2func(funcname);
delete(gcp('nocreate'))
if ~exist('num_workers','var')
  num_workers = 4;
end
parpool(num_workers)
parfor i = 1:numel(JOBs)
  JOB = JOBs{i};
  try
    JOB = func(JOB);
  catch ME
    %{
possible job names:
    fname_epi | myspm_fmriprep12
    filenames{1} | myspm_fmriglm, myspm_glm
    fname | myspm_smooth
    fname_t1w | myspm_seg12, myspm_seg_cat12
    %}
    if isfield(JOB,'fname_epi')
      jobname = JOB.fname_epi;
    elseif isfield(JOB,'filenames')
      jobname = JOB.filenames{1};
    elseif isfield(JOB,'fname')
      jobname = JOB.fname;
    elseif isfield(JOB,'fname_t1w')
      jobname = JOB.fname_t1w;
    else
      jobname = 'unknown';
    end
    warning(ME.identifier, '### ERROR DURING PROCESSING: %s\n', jobname);
    warning(ME.identifier, '### %s', ME.message);
    for iStack=1:numel(ME.stack)
      fprintf('### In %s (line %i):\n', ...
        ME.stack(iStack).file, ME.stack(iStack).line);
    end
    JOB.ME = ME;
    JOB.error_with = JOB.fname_epi;
  end
  JOBs{i} = JOB;
end
delete(gcp('nocreate'))
end