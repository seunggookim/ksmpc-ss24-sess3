function TR_msec = ernsttr(fa_deg, B0_Tesla)
% fa_deg = ernstangle(TR_msec, B0_Tesla)
%
%
% SEE ALSO ERNSTANGLE
% (cc) 2019, sgKIM.

if ~exist('B0_Tesla','var')
  B0_Tesla = 3;
end
if B0_Tesla==3
  T1 = 1000; % REF: Weiskopf's paper...
end
TR_msec = -T1*log(cos(fa_deg/180*pi));
end