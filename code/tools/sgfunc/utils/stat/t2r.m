function r = t2r(t,n)
% r = t2r(t,n)
% t = r * sqrt( (n-2) ./ (1-r.^2) )

% ref: https://en.wikipedia.org/wiki/Pearson_correlation_coefficient

r = t ./ sqrt( n-2+t.^2 );
end