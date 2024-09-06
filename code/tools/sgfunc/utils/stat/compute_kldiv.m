function [kl,p1,p2,edges] = compute_kldiv(vec1, vec2, edges)
% computes KL divergence between PDFs of given vectors
% 
% kl = compute_kldiv(vec1, vec2, nbins)
%
% (cc) sgKIM.


% exclude NaNs:
vec1 = vec1(~isnan(vec1));
vec2 = vec2(~isnan(vec2));

% if exist('nbins','var')
%   [n1,edges] = histcounts(vec1,nbins);
% else
%   [n1,edges] = histcounts(vec1);
% end
if exist('edges','var')
  n1 = histcounts(vec1, edges);
else
  [n1,edges] = histcounts(vec1);
end
n2 = histcounts(vec2, edges);
p1 = n1/sum(n1);
p2 = n2/sum(n2);
kl = sum(p1.*log((p1+eps)./(p2+eps)));
end