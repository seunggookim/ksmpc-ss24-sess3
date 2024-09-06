function [t,p] = r2t(r,n,tail)
% [t,p] = r2t(r,n,tail)
%
% INPUTS
% r : Pearson correlation coefficient(s)
% n : # of samples
% tail : 'both' (default) | 'left' | 'right'
%
% OUTPUTS
% t : t-stat(s)
% p : p-value(s) with the given tail(s)
%
% This is basically a wrapper of PVALPEARSON from CORR.
% (cc) 2019, sgKIM. solleo@gmail.com

if ~exist('tail','var')
  tail = 'b';
end
[p,t] = pvalPearson(tail, r, n);
end

function [p,t] = pvalPearson(tail, rho, n)
% p = pvalPearson(tail, rho, n)
% PVALPEARSON Tail probability for Pearson's linear correlation.
% - extracted from CORR(); Copyright 1984-2018 The MathWorks, Inc.

t = rho.*sqrt((n-2)./(1-rho.^2)); % +/- Inf where rho == 1
switch lower(tail(1))
case 'b' % 'both or 'ne'
    p = 2*tcdf(-abs(t),n-2);
case 'r' % 'right' or 'gt'
    p = tcdf(-t,n-2);
case 'l' % 'left' or 'lt'
    p = tcdf(t,n-2);
end
end
