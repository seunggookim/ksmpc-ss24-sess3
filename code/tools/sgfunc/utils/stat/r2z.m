function [z, p] = r2z(r, n)
% [z, p]=r2z(r, n)
%
% (cc) 2018-2019, sgKIM

% REF: https://en.wikipedia.org/wiki/Fisher_transformation
% z = arctanh(r) Z~N(0,SE=1/sqrt(n-3))

z = atanh(r)*sqrt(n-3);
p = 2*normcdf(-abs(z),0,1); % both tails
end