function Centroids = centroiddist(X)
%CENTROIDDIST finds column-wise centroids.
% Centroids = centroiddist(X)
% X can be either (1) a column vector or (2) a matrix that describes a discrete distribution.


assert(not(isRowVectorWithLength(X)), 'X can be either (1) a column vector or (2) a matrix!')

if min(X(:))<0
  warning('Input X includes negative values. Adjusting the offset to make them nonnegative.')
  X = X - min(X(:));
end

Idx = (1:size(X,1))';
WeightedSum = sum(repmat(Idx,[1, size(X,2)]) .* X, 1);
Centroids = WeightedSum ./ sum(X,1);

end