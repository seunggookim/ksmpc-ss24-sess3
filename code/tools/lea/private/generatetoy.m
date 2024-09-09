function [X, Y] = generatetoy(Job)
% [X, Y] = generatetoy(Job)

X = {};
Y = {};
for iSet = 1:Job.nSets
  X_ = normrnd(0, 1, [Job.nSamples, Job.nFeatures]);
  Y_ = normrnd(0, 1, [Job.nSamples, Job.nResponses]);
  Y_(:,1) = Y_(:,1) + X_(:,1)*Job.EffectSize;
  if Job.TempGaussWin
    X_ = smoothdata(X_, 'gauss', Job.TempGaussWin);
    Y_ = smoothdata(Y_, 'gauss', Job.TempGaussWin);
  end
  X = [X, X_];
  Y = [Y, Y_];
end
end
