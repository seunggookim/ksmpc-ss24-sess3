function Job = setdefaultjob(DefaultJob, Job)
%SETDEFAULTJOB check input fields and set default values
%
% Job = setdefaultjob(DefaultJob, Job)
%
% (cc) 2021, sgKIM.

st = dbstack;
ProcName = st(2).name;

FldNames = fieldnames(DefaultJob);
for iFld = 1:numel(FldNames)
  if ~isfield(Job, FldNames{iFld})
    fprintf('[%s:%s] (DEFAULT) Job.%s = ', ...
      ProcName, datestr(now,31), FldNames{iFld});
    disp(DefaultJob.(FldNames{iFld}))
    Job.(FldNames{iFld}) = DefaultJob.(FldNames{iFld});
  end
end

end
