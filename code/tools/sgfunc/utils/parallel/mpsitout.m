function mpsitout(FuncHandle, Jobs)
%mpsitout(FuncHandle, Jobs)
%
s = functions(FuncHandle);
assert(~isempty(s.file), '"%s" is not pathed.', FuncHandle)
assert(iscell(Jobs), 'Jobs should be a cell array!')

tStart = tic;
FnameSh = [tempname,'_runme.sh'];

% but how can you 


system(['sbatch --wait ',FnameSh])
fprintf('[%s:%s] Done: ', mfilename, datestr(now,31)); toc(tStart)


end