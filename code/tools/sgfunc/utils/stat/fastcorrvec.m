function [r] = fastcorrvec(A,B)
% [r] = fastcorrvec(A,B)
%
% computes a correlation vector r (1xP) between matrices A (NxP) and B (NxP)
%
% fast (~20x than CORR; 1.5x than CORRVEC) but slightly off (err < 1e-14)
%
% (cc) 2023, dr.seunggoo.kim@gmail.com

% REF: https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
%
%                          n*sum(x_i,y_i) - sum(x_i)*sum(y_i)
% r_xy = ------------------------------------------------------------------
%         sqrt(n*sum(x_i^2)-(sum(x_i))^2) * sqrt(n*sum(y_i^2)-(sum(y_i))^2)
%

n = size(A,1);
mA = mean(A,1);
mB = mean(B,1);
Nom = bsxfun(@minus, sum(bsxfun(@times, A, B)),  n*bsxfun(@times, mA, mB));
Denom = sqrt( bsxfun(@times, ...
  bsxfun(@minus, sum(A.^2), n*mA.^2), bsxfun(@minus, sum(B.^2), n*mB.^2) ) );
r = bsxfun(@rdivide, Nom, Denom);
end