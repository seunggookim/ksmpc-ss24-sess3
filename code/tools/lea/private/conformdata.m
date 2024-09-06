function Data = conformdata(X, Y, Job)
% Data = conformdata(X, Y, Job)
% (cc) 2024, seung-goo.kim@ae.mpg.de

assert(numel(X) == numel(Y), "Dimensions of X and Y mismatched")
nSets = numel(X);
Data = [];

for iSet = 1:nSets
  % Find out overlapping time samples and apply the relative time of interest
  TimesX = (0:size(X{iSet},1)-1)'/Job.SamplingRateHz;
  TimesY = (0:size(Y{iSet},1)-1)'/Job.SamplingRateHz;
  MaxTimeSec = min(TimesX(end), TimesY(end));
  TimeMaskX = (Job.RelToiSec(1) <= TimesX) & ...
    (TimesX <= (MaxTimeSec + Job.RelToiSec(2)));
  TimeMaskY = (Job.RelToiSec(1) <= TimesY) & ...
    (TimesY <= (MaxTimeSec + Job.RelToiSec(2)));
  
  % Delay, crop, & standardize stimulus
  X_ = zscore(delayreg(X{iSet}, TimeMaskX, Job.DelaysSmp, false));
  X_ = [ones(size(X_,1),1), X_]; % adding a bias term

  % Crop & standardize response
  Y_ = zscore(Y{iSet}(TimeMaskY,:));

  % Sample times
  T_ = TimesX(TimeMaskX);

  % Contain all in a structure:
  Data = [Data, struct(X=X_, Y=Y_, T=T_)];
end

end
