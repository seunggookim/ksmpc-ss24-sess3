function [stat, cfg] = myft_permt2(dat,group,cov)
% [s, cfg] = myft_permt2(y,group,cov)

nSubj = numel(group);
if size(dat,2) ~= nSubj
  dat = dat';
end
if iscolumn(group)
  group = group';
end
if islogical(group) || all(unique(group)==[0 1])
  group = group+1;
end
if exist('cov','var')
  dat = residual(dat', double(term(cov)+1))';
end

design = group;
cfg = struct('alpha',0.05, 'tail',0, 'correcttail','prob', ...
  'ivar',1, 'statistic','indepsamplesT', 'numrandomization',1e+5, ...
  'resampling','permutation','permute_only_ivar',1);
[stat, cfg] = ft_statistics_montecarlo(cfg, dat, design);
end