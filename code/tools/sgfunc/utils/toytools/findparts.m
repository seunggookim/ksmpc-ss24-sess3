function idx = findparts(n,k)
% idx = findparts(n,k)
%
% For a sequence [1,2,...n], create k sets comprising all elements of the
% sequence.
%
% (cc) 2022, dr.seunggoo.kim@gmail.com

p = floor(n/k);
idx = arrayfun(@(x,y) x:(x+y-1), [1:p:(p*k)], p*ones(1,k), 'uni',0);
idx = [idx (p*k+1):n];
end