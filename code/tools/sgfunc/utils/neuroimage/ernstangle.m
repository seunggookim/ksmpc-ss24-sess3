function fa_deg = ernstangle(TR_msec, B0_Tesla)
% fa_deg = ernstangle(TR_msec, B0_Tesla)
%
%
% SEE ALSO ERNSTTR
% (cc) 2019, sgKIM.

if ~exist('B0_Tesla','var')
  B0_Tesla = 3;
end
if B0_Tesla==3
  T1 = 1000; % REF: Weiskopf's paper...
end
fa_deg = acos(exp(-TR_msec/T1))/pi*180;
end