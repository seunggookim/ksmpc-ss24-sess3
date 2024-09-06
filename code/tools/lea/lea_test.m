function [Job, Mdl, Rnd] = lea_test(Job)
% [Job, Mdl, Rnd] = lea_test(Job)
% (CC0-BY) seung-goo.kim@ae.mpg.de

% set up parameters🍽️:
if not(nargin); Job = []; end
disp(repmat('=',[1 80]))
Job = defaultjob(struct(nSamples=50, nFeatures=2, nResponses=3, nSets=4, ...
  EffectSize=3, SamplingRateHz=1, DelaysSmp=[0, 1], RelToiSec=[2 -2], ...
  LambdaGrid=10.^(-5:0.5:5), nRands=1000, IsFigure=false, ...
  TempGaussWin=10, IsTest=true), Job);
disp(Job)

% generate toy data🧸:
[X, Y] = generatetoy(Job);

% run cv folds🏃‍♀️‍➡️‍:
Mdl = runcv(X, Y, Job);

logthis('Mean Lopt = ')
disp(geomean(Mdl.Lopt, 1))
fprintf('\b')
if Job.IsTest
  assert( geomean(Mdl.Lopt(:,1), 1) < 10, "OVER-OPTIMIZED" )
  assert( min(geomean(Mdl.Lopt(:,2:end), 1)) > 1, "UNDER-OPTIMIZED" )
end

logthis('Mean acc = ')
disp(mean(Mdl.Acc, 1))
fprintf('\b')
if Job.IsTest
  assert( mean(Mdl.Acc(:,1), 1) > 0.5, "FALSE NEGATIVE" )
  assert( max(mean(Mdl.Acc(:,2:end), 1)) < 0.5, "FALSE POSITIVE")
end

% randomization test👾:
Rnd = randtest(X, Y, Mdl, Job);

% create nice plots📊️:
if Job.IsFigure
  plotmdl(X, Y, Mdl, Rnd, Job)
end

if Job.IsTest
  logthis('ALL PASSED!\n')
end

if not(nargout)
  clear Job
end

end
