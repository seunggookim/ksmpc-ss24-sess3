function a = stat_desc_rm(x1, x2, x1name, x2name)
% a = stat_desc_rm(x1, x2, x1name, x2name) -  Repeated-measures boxplot + scatterplot + individual lines
%
% [!!] changes (differences) of x2 with respect ot x1 (x2 - x1)

% desc stats:
a(1,1) = min(x1);
a(1,2) = max(x1);
a(1,3) = nanmean(x1);
a(1,4) = nanstd(x1);

a(2,1) = min(x2);
a(2,2) = max(x2);
a(2,3) = nanmean(x2);
a(2,4) = nanstd(x2);

a(3,1) = min(x2-x1);
a(3,2) = max(x2-x1);
a(3,3) = nanmean(x2-x1);
a(3,4) = nanstd(x2-x1);

if ~exist('x1name','var'), x1name='x1'; end
if ~exist('x2name','var'), x2name='x2'; end
hold on;
for i=1:numel(x1)
  plot([1 2]',[x1(i) x2(i)]','color',[.7 .7 .7], 'linewidth',.5, 'color',[0.3 0.3 0.3 .5])
end
cfg=struct('showdata',4, 'groupcolor',[0 .8 .1; 1 0 0], 'outliers',0, ...
  'maxwhisker',1.5, 'boxwidth',0.3);
[out,s] = cat_plot_boxplot_sg({x1, x2}, cfg);

[~,p,~,st] = ttest(x2-x1);
title(sprintf('%s-%s, T(%d)=%.2f, P=%.3f',x2name,x1name,st.df,st.tstat,p))
set(gca,'xticklabel',{x1name,x2name}, 'xticklabelrotation',0)
end