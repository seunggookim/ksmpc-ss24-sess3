function [DnLog, JobId] = findslurmlogpath()
%FINDSLURMLOGPATH finds a next natural number for the name of SLURM log folder
% [DnLog, JobId] = findslurmlogpath()

DnTemp = fullfile(getenv('HOME'), 'slurm');
if ~isfolder(DnTemp); mkdir(DnTemp); end
Dirs = dir(fullfile(DnTemp, '*'));
Dirs = Dirs([Dirs.isdir]);
JobId = str2double(Dirs(end).name);
if isnan(JobId)
    JobId = 0; 
else
    JobId = JobId + 1;
end
DnLog = fullfile(DnTemp, sprintf('%03i', JobId));
mkdir(DnLog)
JobId = sprintf('%03i', JobId);
end
