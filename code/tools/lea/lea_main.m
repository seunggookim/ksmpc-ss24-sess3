function Job = lea_main(Job)
% Linearized Encoding Analysis (LEA👧)
%
% (cc) 2024-09-02, seung-goo.kim@ae.mpg.de

% house-keeping🧹:
Job = defaultjob(struct(HighPassHz=0.001, SamplingRateHz=1, ...
  DelaysSmp=[0, 1], RelToiSec=[10 -10], LambdaGrid=10.^(-5:0.5:5), ...
  nRands=1000, DnameMdl=tempdir), Job);
[~,~] = mkdir(Job.DnameMdl);
FnameLog = fullfile(Job.DnameMdl, 'lea_main.log');
tStart = tic;
logfile(FnameLog, Job, true)

% preprocess data🍱:
[X, Y] = prepdata(Job);

% run cv folds🏃‍♀️‍➡️‍:
Mdl = runcv(X, Y, Job);
Job.FnameMdl = fullfile(Job.DnameMdl, 'mdl.mat');
save(Job.FnameMdl, 'Mdl', 'Job')

% randomization test👾:
Rnd = randtest(X, Y, Mdl, Job);
Job.FnameRnd = fullfile(Job.DnameMdl, 'rnd.mat');
save(Job.FnameRnd, 'Rnd', 'Job')

% create nice plots📊️:
plotmdl(X, Y, Mdl, Rnd, Job)

% close the log📝:
logfile(FnameLog, Job, false, tStart)
end
