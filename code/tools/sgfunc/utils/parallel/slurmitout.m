function slurmitout(FuncHandle, Jobs, InitCmd, InitSh, Cfg)
%slurmitout(FuncHandle, Jobs, InitCmd, InitSh, Cfg)
%
% e.g.
% InitCmd = 'addpath /my/local/path; myfunc_addpath';
% InitSh  = 'source ~/.bashrc';
%
% FuncHandle (1x1)      MATLAB function handle
% Jobs       {nJobsx1}  Jobs in a cell array
% InitCmd    '1xnChar'  A one-lined MATLAB command to run before run the given function (default='')
% InitSh     '1xnChar'  A one-lined Shell command to run before open MATLAB (default='source ~/.bashrc')
%
% Cfg        [1x1]      SLURM configuration switches
% .nTasks     [1x1]      The number of tasks per job (default=1)
% .Mem_GB     [1x1]      Memory per node in GB (default=8)
% .Partition  '1xnChar'  Partition name in SLURM (default='octopus')
% .IsWait     [T\F]      --wait switch (default=true)
% .Switch     '1xnChar'  Additional SLURM switches
%
% (cc) 2022-2023, dr.seunggoo.kim@gmail.com

s = functions(FuncHandle);
assert(~isempty(s.file), '"%s" is not pathed.', FuncHandle)
assert(iscell(Jobs), 'Jobs should be a cell array!')
tStart = tic;
if ~exist('InitCmd','var')||isempty(InitCmd); InitCmd=''; end
if ~exist('InitSh','var')||isempty(InitSh); InitSh='source ~/.bashrc'; end
if ~exist('Cfg','var')||isempty(Cfg); Cfg = struct(); end
[DnTemp, JobId] = findslurmlogpath();

%% Save Job files
for iJob = 1:numel(Jobs)
  Job = Jobs{iJob};
  FnameJob = [DnTemp,filesep,'Job_',num2str(iJob),'.mat'];
  save(FnameJob,'Job')
end

%% Create runme.sh
FnameSh = fullfile(DnTemp,'runme.sh');
DefaultCfg = struct('Partition','octopus', ...
  'CpuPerTask', 1, 'nTasks', 1, 'Mem_GB', 5, ...
  'IsWait', true, 'Switch', '');
Cfg = defaultcfg(DefaultCfg, Cfg, mfilename);
fprintf('\n')

fid = fopen(FnameSh,'w');
fprintf(fid, '#!/bin/bash\n');
fprintf(fid, '#SBATCH --job-name=%s\n', JobId);
FuncPath = fileparts(s.file);
fprintf(fid, '#SBATCH --chdir=%s\n', FuncPath);
fprintf(fid, '#SBATCH --output=%s/%s-%%A_%%a.out\n', DnTemp, s.function);
fprintf(fid, '#SBATCH --array=1-%i\n', numel(Jobs));
fprintf(fid, '#SBATCH --partition=%s\n', Cfg.Partition);
fprintf(fid, '#SBATCH --ntasks=%i\n', Cfg.nTasks);
fprintf(fid, '#SBATCH --cpus-per-task=%i\n', Cfg.CpuPerTask);
fprintf(fid, '#SBATCH --mem=%iG\n', Cfg.Mem_GB);
fprintf(fid, '#SBATCH -t 14-0:00 %s\n', Cfg.Switch); % time limit: 14 days
fprintf(fid, '%s\n', InitSh);
fprintf(fid, ...
  ['matlab -nodisplay -r "%s; ',...          % open MATLAB
  'load %s_${SLURM_ARRAY_TASK_ID}.mat; ',... % load the i-th job file
  '%s(Job); exit;"\n'], ...                  % run that job file
  InitCmd, [DnTemp,filesep,'Job'], s.function);
fclose(fid);

%% Submit it
logthis('Running %s: %i jobs via SLURM\n', ...
  func2str(FuncHandle), numel(Jobs))
if Cfg.IsWait
  system(['sbatch --wait ',FnameSh]);
  logthis('Done: '); toc(tStart)
else
  system(['sbatch ',FnameSh]);
end
logthis('Check out SLURM outputs: "%s"\n\n', DnTemp);

end
