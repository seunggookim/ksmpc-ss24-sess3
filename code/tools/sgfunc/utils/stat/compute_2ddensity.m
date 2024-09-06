function F = compute_2ddensity(Xi, Y, varargin)

[X1, X2] = ndgrid(Xi, Xi);
F = mvksdensity(Y, [X1(:), X2(:)], varargin{:});
F = reshape(F, size(X1))';
end