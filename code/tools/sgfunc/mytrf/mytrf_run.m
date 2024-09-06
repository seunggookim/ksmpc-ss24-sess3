function job = mytrf_run(job)

if ~isa(job.prepfunc,'function_handle')
  error('job.prepfunc SHOULD BE a function handle to prepare STIM and RESP')
end

if ~isfield(job,'nthreads'), job.nthreads = 16; end
if isfield(job,'conds')
  suffix = ['cnd-',cell2mat(job.conds)];
elseif isfield(job,'randstim')
  suffix = ['rnd-',job.randstim];
else
  suffix = ['cnd-all'];
end

if isfield(job,'exceptsub')
  job.subID = setdiff(job.subID, job.exceptsub);
end

% RCA:
if isfield(job,'rca')
  myrca_run(job)
end

% set TRF runs:
jobs = cell(numel(job.subID),1);
for i = 1:numel(job.subID)
  jobs{i} = job;
  jobs{i}.subid = job.subID{i};
  if isfield(job,'exceptruns')
    jobs{i}.exceptruns = job.exceptruns{i};
  end
end

% for test run:
if isfield(job,'testrun') && job.testrun
  warning('**TEST RUN**: RUNNING ONLY THE FIRST JOB')
  testjob = runthiscv(jobs{1});
  return
end

% otherwise batch run:
delete(gcp('nocreate'))
parpool(job.nthreads)
parfor i = 1:numel(job.subID)
  jobs{i} = runthiscv(jobs{i});
end
delete(gcp('nocreate'))

% Create figure:
d = dir(strrep(jobs{1}.fn_cv, job.subID{1},'*'));
if isfield(job,'exceptsub')
  d(contains({d.folder}, job.exceptsub)) = [];
  suffix = [suffix,'_exceptsub'];
end
job = [];
job.fnames = strcat({d.folder}', filesep, {d.name}');
% just enclose different analyses together:
job.dn_out = fileparts(...
  strrep(strrep(jobs{1}.fn_cv, jobs{1}.subid, 'figure'),['_',suffix],''));
job.dn_out = fileparts(strrep(jobs{1}.fn_cv, jobs{1}.subid, 'figure'));
job.suffix = suffix;
mytrf_plotgrouptrf(job)

end

function job = runthiscv(job)
% stim/data prep
[STIM,RESP,job] = feval(job.prepfunc, job);

% model opt+eval (nested CV)
[~,job] = mytrf_cv(STIM,RESP,job);
end
