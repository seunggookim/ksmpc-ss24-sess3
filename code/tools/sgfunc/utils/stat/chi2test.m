function [chi2 ,pval, df, p_h0, m_obs, m_exp]= chi2test (classID, groupID)
%CHI2TEST tests the equoprobable K classes across N groups.
% [chi2, pval, df, p_h0, m_obs, m_exp] = chi2test (label, groupID)
%
% [EXAMPLES]
% gender difference (1=female, 2=male) across 3 groups
% [chi2, pval, df] = chi2test([1 2 1 2 1 2]',[1 1 2 2 3 3]')
%
% language difference (1=eng, 2=chn, 3=kor) across 2 groups
% [chi2, pval, df] = chi2test([1 2 3 1 2 3]',[1 1 1 2 2 2]')
%
% SEE ALSO: chi2cdf
%
% (cc) sgKIM, 2012, 2017

%{
REF: https://web.stanford.edu/class/psych10/schedule/P10_W7L1
General approach:
- Generate a null hypothesis (H0) about population proportions.
- Compute frequencies that we would expect if H0 was true.
- Summarize the discrepancy between these with a single statistic (Χ2)
- Use the distribution of all of the statistics that we could have observed
if the null hypothesis was true to determine whether this statistic would
be unlikely if the null hypothesis was true.
%}

% CHECK inputs:
assert(isvector(classID), 'CLASSID should be a vector!')
assert(isvector(groupID), 'GROUPID should be a vector!')

% MAKE binary labels for classes
classes = unique(classID);
nclasses = numel(classes);
binlabels = false(numel(classID), nclasses);
for iclass = 1:nclasses
  binlabels(classID==classes(iclass),iclass) = true;
end

% H0: equiprobable across groups
p_h0 = mean(binlabels);

% m: # of samples in the contingency table
groups = unique(groupID);
ngroups = numel(groups);
m_obs = zeros(nclasses,ngroups);
m_exp = zeros(nclasses,ngroups);
for iclass = 1:nclasses
  for igroup = 1:ngroups
    ind_group = groupID==groups(igroup);
    m_obs(iclass, igroup) = sum(binlabels(ind_group, iclass));
    m_exp(iclass, igroup) = sum(ind_group) * p_h0(iclass);
  end
end

% COMPUTE stats:
chi2 = sum( (m_obs(:)-m_exp(:)).^2 ./ m_exp(:) );
df = (ngroups-1) * (nclasses-1);
pval = 1 - chi2cdf(chi2, df);


end
