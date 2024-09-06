function slurmitoutinbatch(FuncHandle, Jobs, varargin)
%slurmitinbatch(FuncHandle, Jobs, InitCmd, InitSh, Cfg)
%
% SEE ALSO: SLURMITOUT
% (cc) 2022-2023, dr.seunggoo.kim@gmail.com

BatchSize = varargin{3}.BatchSize;
nBatch = ceil(numel(Jobs)/BatchSize);
iBatch = 1;
while not(isempty(Jobs))
  % get the first [Cfg.BatchSize] (or less) jobs
  Idx = 1:min(BatchSize, numel(Jobs)); 
  JobsBatch = Jobs(Idx); % have the subset in a new array
  Jobs(Idx) = []; % remove it from the array
  logthis('Running %i/%i batch...\n', iBatch, nBatch)
  slurmitout(FuncHandle, JobsBatch, varargin{:})
  iBatch = iBatch + 1;
end

end
