Jobs = {};
for i = 1:1000
  Job = struct('ArrayClass','double','nSub',30, 'nRnd',2000, 'iDraw',i);
  Job.DnSim = ['~/usrfolder/clusterstatminp_test/',...
    sprintf('n%i_k%i', Job.nSub, Job.nRnd)];
  Jobs = [Jobs Job];
end
InitCmd = ['addpath /home/seung-goo.kim/git/NCML-code/matlab/sgfunc;' ...
  'addpath_sgfunc;'];
slurmitout(@clusterstatminp_test, Jobs, InitCmd)
%[slurmitout:2022-01-04 14:29:14] Done: Elapsed time is 234.465293 seconds.

%%
Files = dir(['~/usrfolder/clusterstatminp_test/n30_k2000/*.mat']);
Fnames = strcat({Files.folder},'/',{Files.name});
minp = []; posminp = []; maxt = [];
for i = 1:numel(Fnames)
  Data_ = load(Fnames{i}, 'stat','StatObs');
  minp(i) = min(Data_.stat.prob);
  posminp(i) = min(Data_.stat.posobsminp);
  maxt(i) = max(Data_.StatObs(:));
end
figure;
histogram(minp)
% okay... kind of flat?
