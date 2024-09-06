function Job = defaultjob(DefaultJob, Job, ProcName)
%defaultjob check input fields and set them to default values
%
% Job = defaultjob(DefaultJob, Job)
%
% (cc) 2021, sgKIM.

if not(exist('ProcName','var')); ProcName=[]; end
FldNames = fieldnames(DefaultJob);
for iFld = 1:numel(FldNames)
  if ~isfield(Job, FldNames{iFld})
    if not(isempty(ProcName))
      fprintf('[%s] (DEFAULT) Job.%s = ', ProcName, FldNames{iFld})
      disp(DefaultJob.(FldNames{iFld})); fprintf('\b')
    end
    Job.(FldNames{iFld}) = DefaultJob.(FldNames{iFld});
  end
end

end
