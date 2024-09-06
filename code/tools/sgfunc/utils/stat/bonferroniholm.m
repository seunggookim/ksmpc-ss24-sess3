function alpha_adj = bonferroniholm(pvals, alpha0)
% Computes adjusted alpha with Bonferroni-Holm correction
%
% USAGE:
% alpha_adj = bonferroniholm(pvals, alpha)
% 
% REF: Cramer et al., 2016, Psych Bull Rev
%
% (cc) 2020, sgKIM, solleo@gmail.com

if ~exist('alpha','var')
  alpha0 = 0.05;
end

% Adjusted alpha = alpha / reverse-order of p-value
[~,idx] = sort(pvals);
[~,order] = ismember(pvals, pvals(idx));
revorder = numel(pvals)-order+1;
alpha_adj = alpha0./revorder;

% Let all H0 accepted after the first acceptance:
k = find(pvals(idx) > alpha_adj(idx), 1, 'first');
if ~isempty(k)
  alpha_adj(idx(k:end)) = -inf;
end

end