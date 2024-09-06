function d = compute_dprime(isHit,isFalsealarm, method)
% d = compute_dprime(isHit,isFalsealarm, method)
%
% extreme values are corrected based on Macmillan & Kalplan (1985)
%
% REF-1: https://en.wikipedia.org/wiki/Sensitivity_index
% REF-2: https://stats.stackexchange.com/questions/134779/d-prime-with-100-hit-rate-probability-and-0-false-alarm-probability
% METHOD-1(*): Adjust only the extreme values by replacing rates of 0 with 0.5/n and rates of 1 with (n-0.5)/n where n is the number of signal or noise trials (Macmillan & Kaplan, 1985, Detection theory analysis of group data: estimating sensitivity from average hit and false-alarm rates, Psychol Bull.)
% MEHOTD-2: add 0.5 to both the number of hits and the number of false alarms, and add 1 to both the number of signal trials and the number of noise trials; dubbed the loglinear approach (Hautus, 1995)
%
% (cc) 2019-02-15, sgKIM, solleo@gmail.com.

if ~exist('method','var'), method = 1; end

switch method % "Correction" of extreme values (0 or 1) to avoid �Inf
  case 1 % METHOD-1 Macmillan & Kaplan
    hitrate = mean(isHit);
    n1 = numel(isHit);
    farate = mean(isFalsealarm);
    n2 = numel(isFalsealarm);
    hitrate = hitrate +(0.5/n1)*(hitrate==0) -(0.5/n1)*(hitrate==1);
    farate = farate +(0.5/n2)*(farate==0) -(0.5/n2)*(farate==1);
  case 2 % METHOD-2 Hautus
    hitrate = (sum(isHit)+0.5)/(numel(isHit)+1);
    farate= (sum(isFalsealarm)+0.5)/(numel(isFalsealarm)+1);
  otherwise
    error('unknown method?')
end

d = icdf('norm',hitrate,0,1) - icdf('norm',farate,0,1);
end