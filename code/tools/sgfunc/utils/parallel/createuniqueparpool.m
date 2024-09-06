function p = createuniqueparpool(nWorkers)
% randomize job storage name (so that I can run these in parallel)
% Pool = createuniqueparpool(nWorkers)
%
% Solution provided by the MATLAB support
% (cc) 2023, dr.seunggoo.kim@gmail.com

c = parcluster('local');
c.NumWorkers = nWorkers;
Pid = java.lang.management.ManagementFactory.getRuntimeMXBean.getName.char;
JobPath = fullfile(c.JobStorageLocation, Pid, getenv('SLURM_ARRAY_TASK_ID'));
[IsDone, ErrMsg, ErrId] = mkdir(JobPath);
if not(IsDone)
  error(ErrId, 'A new folder "%s" not created: %s', JobPath, ErrMsg)
end
c.JobStorageLocation = JobPath;
delete(gcp('nocreate'));
p = parpool(c);
end
