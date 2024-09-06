function norm = l2norm(x)
% norm = sqrt ( sum ( x.^2, 2) )
% norm = sqrt ( sum ( x.^2, 2) );
warning('Don''t re-invent the wheel: use VECNORM(X,2,2) for l2-norm along the columns (leaving rows intact)!');
x = vecnorm(x,2,2);
end
