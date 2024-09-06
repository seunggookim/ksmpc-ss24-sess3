function [X, Y] = generatetoy(Job)
% [X, Y] = generatetoy(Job)

X = {};
Y = {};
for iSet = 1:Job.nSets
  X_ = normrnd(0, 1, [Job.nSamples, Job.nFeatures]);
  Y_ = normrnd(0, 1, [Job.nSamples, Job.nResponses]);
  if Job.TempGaussWin
    X_ = smoothdata(X_, 'gauss', Job.TempGaussWin);
    Y_ = smoothdata(Y_, 'gauss', Job.TempGaussWin);
  end
  X = [X, X_];
  Y = [Y, Y_];
  Y{end}(:,1) = Y{end}(:,1) + X{end}(:,1)*Job.EffectSize;
end
end
