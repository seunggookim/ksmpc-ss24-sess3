function [mR,R] = compute_isc(Y)
% computes inter-subject correlation
%
% [mR,R] = compute_isc(Y)
% Y: #times x #voxels x #subjects
% mR: 1 x #voxels
% R:  1 x #voxels x #subjects

[~, nvoxs, nsubs] = size(Y);
R = nan([1 nvoxs nsubs], class(Y));
fprintf('Computing ISC')
for isub = 1:nsubs
  mY = mean(Y(:,:,setdiff(1:nsubs, isub)),3);
  R(:,:,isub) = corrvec(mY, Y(:,:,isub));
  fprintf('.')
end
fprintf('\n')
mR = mean(R,3);
end



function [r,p] = corrvec(A,B,tail)
% [r,p] = corrvec(A,B,[tail])
%
% computes a correlation vector r (1xP) between matrices A (NxP) and B (NxP)
%
% (cc) 2019, sgKIM. solleo@gmail.com
%
% - THIS CAN BE A MEMORY SURGE...

% REF: https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
%                   sum(x_i - x_bar)(y_i - y_bar)
% r_xy = -------------------------------------------------
%        sqrt(sum(x_i-x_bar)^2) . sqrt(sum(x_i-x_bar)^2) )

if ~exist('tail','var'), tail='b'; end
A = A - mean(A,1);
B = B - mean(B,1);
r = sum(A.*B,1);
r = r ./ ( sqrt(sum(A.^2,1)) .* sqrt(sum(B.^2,1)) );
p = pvalPearson(tail, r, size(A,1)); % using T-distribution
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


function r2 = r2vec(Y,Yhat)
% r2 = r2vec(y,yhat)
%
% computes R2 for each column vector of Y (NxP) given Yhat (NxP):
% Y = [y1, y2, ..., yp]
% y1 = [y11, y12, ..., y1n]'
%
%             sum( (yi - yi_hat)^2 )
% R2_i = 1 - ------------------------
%             sum( (yi - yi_bar)^2 )
%
% where yi_bar is a mean of yi, yi_hat is a prediction for yi.
%
% (cc) 2019, sgKIM. solleo@gmail.com

if ~(all(size(Y)==size(Yhat)))
  error('Y and Yhat should be in the same dimension!')
end
r2 = 1 - sum((Y-Yhat).^2)./sum((Y-ones(size(Y)).*mean(Yhat,1)).^2);
end
