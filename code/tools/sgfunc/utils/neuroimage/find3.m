function [coords, vals] = find3(X,Y)
% *** this is just a wrapper of find3.m by Simon Robinson ***
% http://www.mathworks.com/matlabcentral/fileexchange/30853-field-mapping-toolbox/content//find3.m
%
% [coords] = find3(X)
% returns coordinates of non-zero values in 3-dimensional matrix X
%
% [coords, vals] = find3(X,Y)
% returns coordinates of non-zero values in 3-dimensional matrix X 
% and extracts values in 3-dimensional matrix Y (identical dimenson as X)
%
% Modified by sgKIM.


idx = find(X);
[i,j,k] = ind2sub(size(X), idx);
coords=[i j k];
if nargin > 1
vals = Y(idx);
else
vals = [];
end

end
