function R = rescalecol(A, varargin)
%RESCALECOL Rescale data per each column
% R = rescalecol(A, varargin)
%
% SEE ALSO: RESCALE
%
% (cc) 2020, sgKIM.

R = nan(size(A));
for icol = 1:size(A,2)
  R(:,icol) = rescale(A(:,icol), varargin{:});
end
end