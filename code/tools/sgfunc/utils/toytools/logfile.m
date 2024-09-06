function logfile(FnameLog, Job, IsStart, tStart)
%LOGFILE creates or closes a log file.
%
% [] = logfile(FnLog, Job, IsStart, tStart)
%
% (cc) 2021, sgKIM.

[~, ProcName,~] = fileparts(FnameLog);
DateStrNow = char(datetime('now'),'yyyy-MM-dd''T''HH:mm:ss');
% DateStrNow = char(datetime('now'),'yyyy-MM-dd HH:mm:ss');

if exist('IsStart','var') && IsStart
  diary(FnameLog)
  disp(repmat('=',[1 72]))
  fprintf('[%s|%s] START\n', ProcName, DateStrNow)
  disp(Job)
  save(strrep(FnameLog,'.log','.mat'), 'Job')
else
  if exist('tStart','var')
    fprintf('[%s|%s] END: took %.3f sec\n', ProcName, DateStrNow, toc(tStart))
  else
    fprintf('[%s|%s] END\n', ProcName, DateStrNow)
  end
  diary off
  save(strrep(FnameLog,'.log','.mat'), 'Job')
end
end
