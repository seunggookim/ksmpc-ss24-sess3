function validatefilenames(Job, FldName)
%validatefilenames(Job, FldName)
%
% (cc) 2021, sgKIM.

for iJob = 1:numel(Job)
  Filenames = Job(iJob).(FldName);
  if iscell(Filenames)
    cellfun(@(x) assert(isfile(x), 'FILE NOT FOUND: "%s"',x), Filenames, 'uni',false)

  elseif isstring(Filenames)
    arrayfun(@(x) assert(isfile(x), 'FILE NOT FOUND: "%s"',x), Filenames, 'uni',false)
    
  elseif ischar(Filenames)
    assert(isfile(Filenames), 'FILE NOT FOUND: "%s"',Filenames)
  end
end
end
