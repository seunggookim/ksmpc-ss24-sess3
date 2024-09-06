function slurm1(ShCmd, Mem_GB, Partition)
%slurm1 submit a single job
% slurm1(ShCmd, InitSh, Mem_GB, Partition)
%
% e.g.
% ShCmd = ['source ~/.bashrc; conda activate myshit;',...
% 'python3 /my/python/script /my/argument1 /my/argument2'];
% 
% NOTE1: use /ABSOULTE/PATH for everything!
% NOTE2: SLURM files are saved in ${HOME}/slurm/0xx
%
% (cc) 2022, Seung-Goo KIM


tStart = tic;
DnTemp = fullfile(getenv('HOME'), 'slurm');
if ~isfolder(DnTemp); mkdir(DnTemp); end
Dirs = dir(fullfile(getenv('HOME'), 'slurm', '*'));
JobId = str2double(Dirs(end).name);
if isnan(JobId)
  JobId = 0; 
else
  JobId = JobId + 1;
end
JobId = sprintf('%03i', JobId);
DnTemp = fullfile(getenv('HOME'), 'slurm', JobId);
mkdir(DnTemp)

%% Create runme.sh
FnameSh = fullfile(DnTemp,'runme.sh');
if ~exist('nTasks','var'), nTasks = 1; end
if ~exist('MEM_GB','var'), Mem_GB = 5; end
if ~exist('Partition','var'), Partition='octopus'; end

fid = fopen(FnameSh,'w');
fprintf(fid, '#!/bin/bash\n');
fprintf(fid, '#SBATCH --job-name=%s\n', JobId);
fprintf(fid, '#SBATCH --output=%s/%i.out\n', DnTemp, JobId);
fprintf(fid, '#SBATCH --partition=%s\n', Partition);
fprintf(fid, '#SBATCH --ntasks=1\n');
fprintf(fid, '#SBATCH --cpus-per-task=1\n');
fprintf(fid, '#SBATCH --mem=%iG\n',Mem_GB);
fprintf(fid, '%s\n', ShCmd);
fclose(fid);

%% Submit it
fprintf('[%s:%s] Submitting: ', mfilename, datestr(now,31))
system(['sbatch --wait ',FnameSh]);
fprintf('[%s:%s] Done: ', mfilename, datestr(now,31)); toc(tStart)
fprintf('[%s:%s] Check out SLURM outputs: "%s"\n\n', ...
  mfilename, datestr(now,31), DnTemp);

end
