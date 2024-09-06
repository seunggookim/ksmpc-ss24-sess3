function slurmitoutuntildone(FnamesToDo, varargin)
%slurmitoutuntildone(FnamesToDo, FuncHandle, Jobs, InitCmd, InitSh, Cfg)
%
% SEE ALSO: SLURMITOUT
% (cc) 2022-2023, dr.seunggoo.kim@gmail.com


iCount = 1;
while not(all(cellfun(@isfile, FnamesToDo))) && (iCount<5)
%   tic
  slurmitoutinbatch(varargin{:})
%   Time = toc;
%   logthis('[iCount=%i] Took %f seconds.\n', iCount, Time)
  logthis('[iCount=%i] DONE: %.2f%%\n', ...
    iCount, 100*mean(cellfun(@isfile, FnamesToDo)))
  iCount = iCount + 1;
end
