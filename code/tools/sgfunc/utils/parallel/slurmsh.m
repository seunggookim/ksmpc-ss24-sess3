function slurmsh(ShCmd, Cfg)
%slurmsh submit a single-line (piped, semicoloned, whatever) Bash command
%
% [SYNTAX]
% slurmsh(ShCmd, Cfg)
%
% [INPUT]
% ShCmd       '1xnChar'  A one-lined Shell command to run
% Cfg         [1x1]      SLURM configuration switches
% .nTasks     [1x1]      The number of tasks per job (default=1)
% .Mem_GB     [1x1]      Memory per node in GB (default=8)
% .Partition  '1xnChar'  Partition name in SLURM (default='octopus')
% .IsWait     [T\F]      --wait switch (default=true)
% .Switch     '1xnChar'  Additional SLURM switches
%
% [EXAMPLES]
% e.g.
% ShCmd = ['source ~/.bashrc; conda activate myshit;',...
% 'python3 /my/python/script /my/argument1 /my/argument2'];
% 
% NOTE1: use /ABSOULTE/PATH for everything!
% NOTE2: SLURM files are saved in ${HOME}/slurm/0xx
%
% (cc) 2022-2023, Seung-Goo KIM

assert(exist('ShCmd','var'), 'ShCmd [shell command] is not defined!')
tStart = tic;
[DnTemp, JobId] = findslurmlogpath();
if ~exist('Cfg','var')||isempty(Cfg); Cfg = struct(); end

%% Create runme.sh
FnameSh = fullfile(DnTemp,'runme.sh');
DefaultCfg = struct('Partition','octopus', 'CpuPerTask', 1, 'nTasks', 1, 'Mem_GB', 5, 'IsWait', false, 'Switch', '');
Cfg = defaultcfg(DefaultCfg, Cfg, mfilename);

fid = fopen(FnameSh,'w');
fprintf(fid, '#!/bin/bash\n');
fprintf(fid, '#SBATCH --job-name=%s\n', JobId);
fprintf(fid, '#SBATCH --output=%s/%s.out\n', DnTemp, JobId);
fprintf(fid, '#SBATCH --partition=%s\n', Cfg.Partition);
fprintf(fid, '#SBATCH --ntasks=%i\n', Cfg.nTasks);
fprintf(fid, '#SBATCH --cpus-per-task=%i\n', Cfg.CpuPerTask);
fprintf(fid, '#SBATCH --mem=%iG\n', Cfg.Mem_GB);
fprintf(fid, '#SBATCH -t 14-0:00 %s\n', Cfg.Switch);
fprintf(fid, '%s\n', ShCmd);     % run that shell command
fclose(fid);


%% Submit it
logthis('Running a single job via SLURM\n')
if Cfg.IsWait
  system(['sbatch --wait ',FnameSh]);
  logthis('Done: '); toc(tStart)
else
  system(['sbatch ',FnameSh]);
end
logthis('Check out SLURM outputs: "%s"\n\n', DnTemp);

end
