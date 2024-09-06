function FnLog = matlabnohup(ScriptToRun, PathLog)
% A script to run a MATLAB script via nohup ("no hangups")
%
% [SYNTAX]
% [FnLog] = matlabnohup(ScriptToRun, [PathScript], [PathLog])
%
% [USAGE]
% When you have a script "myScript.m" in a directory in MATLAB's path:
% >> matlabnohup('myScript')
%
% When you have a script "myScript.m" in "/somewhere/else/" that is not 
% in MATLAB's PATH:
% >> matlabnohup('/somewhere/else/myScript.m')
%
% - By default, job logs will be saved in "${HOME}/matlabnohup/" 
%   (e.g., if you use Mac OS and your account is "foo", it will be 
%   "/Users/foo/matlabnohup"). 
% - If you want save the logs in "/whereever/":
%   >> matlabnohup('myScript', '/whereever/')
%
% (cc) 2021, sgKIM. https://github.com/solleo/matlabnohup

if ispc
  error('no nohup in Windows?')
end

%% Set path_log
if ~exist('path_log','var')
  user_home = getenv('HOME');
  PathLog = fullfile(user_home,'matlabnohup');
end
if ~isfolder(PathLog), mkdir(PathLog); end

%% Set jobid
JobId = finduniquejobnumber(PathLog);
PathJob = fullfile(PathLog, JobId);
mkdir(PathJob)

%% Set runme.sh
[PathScript, NameScript,~] = myfileparts(ScriptToRun);
FnSh = fullfile(PathJob,'runme.sh');
Fid = fopen(FnSh,'w');
fprintf(Fid,'#!/bin/bash\nmatlab -nodisplay -r "cd %s; %s; exit"\n', ...
  PathScript, NameScript);
fclose(Fid);

%% Now run:
FnLog = fullfile(PathJob,'run.log');
Fid = fopen(FnLog,'w');
[~,Pid] = system(['nohup bash ',FnSh,' >>',FnLog,' 2>>',FnLog,'& echo $!']);
Pid(end) = []; % remove a linebreak
[~,Hostname] = system('hostname');
Hostname(end) = []; % remove a linebreak
msg = sprintf('[%s:%s] JOB created: [@%s]PID=%s: Logfile=%s\n', mfilename, datestr(now,31), Hostname, Pid, FnLog);
fprintf(Fid,msg);
fclose(Fid);
fprintf(msg);
if ~nargout
  clear FnLog
end
end


function JobId = finduniquejobnumber(pathLog)
% get current jobIDs in the folder:
Files = dir(fullfile(pathLog,'*'));
Files = Files([Files.isdir]);
Files(ismember({Files.name},{'.','..'})) = [];
if ~isempty(Files)
  CurrentIds = cell2mat(cellfun(@str2double, {Files.name}, 'uni',0));
else
  CurrentIds = 0;
end
CurrentIds(isnan(CurrentIds)) = []; % ignore not-a-number
JobId = sprintf('%06i', max(CurrentIds) + 1); % just one bigger
end
