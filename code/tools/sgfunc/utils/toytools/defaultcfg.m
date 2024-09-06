function Cfg = defaultcfg(DefaultCfg, Cfg, ProcName)
%defaultcfg check input fields and set them to default values
%
% Cfg = defaultjob(DefaultCfg, Cfg, ProcName)
%
% (cc) 2021, sgKIM.

FldNames = fieldnames(DefaultCfg);
for iFld = 1:numel(FldNames)
  if ~isfield(Cfg, FldNames{iFld})
    fprintf('[%s] (DEFAULT) Cfg.%s = ', ProcName, FldNames{iFld})
    if isempty(DefaultCfg.(FldNames{iFld}))
      fprintf('\n')
    else
      disp(DefaultCfg.(FldNames{iFld}))
    end
    Cfg.(FldNames{iFld}) = DefaultCfg.(FldNames{iFld});
  end
end

end
