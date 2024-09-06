function [idx,pref] = best(X,dim)
% Finds indices of X that gives maximal value where all values are definite
% and only a single maximum exists (otherwise NaN will be returned).
%
% USAGE
% [idx,pref] = best(X,[dim])
%
% INPUTS
% X   [N_1 x N_2...] a numeric input (vector/matrix/N-D arrays)
% dim [1 x 1]        (optional) a dimension to compute (default = last dim)
%
% OUTPUTS
% idx  [N_dim x 1]
% pref [N_dim x 1] Z-score of the maximum
% (cc) 2020, sgkIM, solleo@gmail.com

if ~exist('dim','var')
  dim = ndims(X);
end

[maxX,idx] = max(X,[],dim);
mask = isnan(sum(X,dim)) | isinf(sum(X,dim)) | (sum(X == maxX,dim) ~= 1);
idx(mask) = nan;
pref = (maxX - mean(X,dim))./std(X,[],dim);
pref(mask) = nan;

end