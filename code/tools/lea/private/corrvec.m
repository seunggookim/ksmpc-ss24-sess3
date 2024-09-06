function [r,p] = corrvec(A,B,tail)
% [r,p] = corrvec(A,B,[tail])
%
% computes a correlation vector r (1xP) between matrices A (NxP) and B (NxP)
%
% (cc) 2019, sgKIM. solleo@gmail.com

% REF: https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
%                   sum(x_i - x_bar)(y_i - y_bar)
% r_xy = -------------------------------------------------
%        sqrt(sum(x_i-x_bar)^2) . sqrt(sum(x_i-x_bar)^2) )

if ~exist('tail','var'), tail='b'; end
deA = A - mean(A,1);
deB = B - mean(B,1);
r = sum(deA.*deB,1) ./ ( sqrt(sum(deA.^2,1)) .* sqrt(sum(deB.^2,1)) );
p = pvalPearson(tail, r, size(A,1)); % using T-distribution
% Z-distribution (Fisher transform) could be also used...
end

function p = pvalPearson(tail, rho, n)
% p = pvalPearson(tail, rho, n)
% PVALPEARSON Tail probability for Pearson's linear correlation.
% - extracted from CORR(); Copyright 1984-2018 The MathWorks, Inc.

t = rho.*sqrt((n-2)./(1-rho.^2)); % +/- Inf where rho == 1
switch tail
case 'b' % 'both or 'ne'
    p = 2*tcdf(-abs(t),n-2);
case 'r' % 'right' or 'gt'
    p = tcdf(-t,n-2);
case 'l' % 'left' or 'lt'
    p = tcdf(t,n-2);
end
end

%{
TESTING CODE:
A = rand(500,100000);
B = rand(500,100000);
clear r p
tic
for i = 1:size(A,2)
  [r(i),p(i)] = corr(A(:,i), B(:,i));
end
toc
tic
[r1,p1] = corrvec(A,B);
toc
mean(r-r1)

Elapsed time is 7.007754 seconds.
Elapsed time is 0.264216 seconds.
ans =
   1.3558e-19
%} 
