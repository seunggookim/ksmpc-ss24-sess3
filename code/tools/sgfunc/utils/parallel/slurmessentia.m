function Job = slurmessentia(Job)

% Check inputs
FieldsNeeded = {'DnPy' 'DnPy' 'FnMdl' 'FnAudio' 'FnOut'};
for i = 1:numel(FieldsNeeded)
  assert(isfield(Job,FieldsNeeded{i}), 'Job."%s" NOT DEFINED!', ...
    FieldsNeeded{i});
end

% Prep a folder to save
DnOut = fileparts(Job.FnOut);
if ~isfolder(DnOut)
  mkdir(DnOut)
end

%% Create runme.sh
tStart = tic;
[DnTemp, JobId] = findslurmlogpath();

FnameSh = fullfile(DnTemp,'runme.sh');
fid = fopen(FnameSh,'w');
fprintf(fid, '#!/bin/bash\n');
fprintf(fid, '#SBATCH --job-name=%s\n', JobId);
fprintf(fid, '#SBATCH --chdir=%s\n', Job.DnPy);
fprintf(fid, '#SBATCH --output=%s/%s.out\n', DnTemp, 'es');
fprintf(fid, '#SBATCH --partition=%s\n', 'gpu');
fprintf(fid, '#SBATCH --ntasks=%i\n',1);
fprintf(fid, '#SBATCH --cpus-per-task=1\n');
fprintf(fid, '#SBATCH --mem=%iG\n',8);
fprintf(fid, 'source ~/.bashrc\n'); % well... it's working...
fprintf(fid, 'conda activate essentia\n');
fprintf(fid, 'python3 %s %s %s %s\n', ...
  Job.FnPy, Job.FnMdl, Job.FnAudio, Job.FnOut);
fclose(fid);

%% Submit it
fprintf('[%s:%s] Submitting: ', mfilename, datestr(now,31))
system(['sbatch --wait ',FnameSh]);
fprintf('[%s:%s] Done: ', mfilename, datestr(now,31)); toc(tStart)
fprintf('[%s:%s] Check out SLURM outputs: "%s"\n\n', ...
  mfilename, datestr(now,31), DnTemp);

%% Check if the output is created
assert(isfile(Job.FnOut), 'File "%s" NOT CREATED!', Job.FnOut)

end
