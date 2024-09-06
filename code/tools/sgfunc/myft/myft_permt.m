function [stat, cfg] = myft_permt(dat,reg,cov)
% [s, cfg] = myft_permt(y,reg,cov)

nSubj = numel(reg);
if size(dat,2) ~= nSubj
  dat = dat';
end
if iscolumn(reg)
  reg = reg';
end
if exist('cov','var')
  dat = residual(dat', double(term(cov)+1))';
end
design = reg;
cfg = struct('alpha',0.05, 'tail',0, 'correcttail','prob', ...
  'ivar',1, 'statistic','indepsamplesregrT', 'numrandomization',1e+5, ...
  'resampling','permutation','permute_only_ivar',1);
[stat, cfg] = ft_statistics_montecarlo(cfg, dat, design);
end