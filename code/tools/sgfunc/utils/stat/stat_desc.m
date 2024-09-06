function a = stat_desc(x,g, groupname, xname)
% a = stat_desc(x,g, groupname, xname)

if nargin == 1
  a(1,1) = min(x);
  a(1,2) = max(x);
  a(1,3) = nanmean(x);
  a(1,4) = nanstd(x);
  cfg = struct('showdata',3, 'groupcolor',[0 .8 .1], 'outliers',0,'maxwhisker',1.5);
  cat_plot_boxplot_sg(x, cfg);
  set(gca,'xtick',1,'xticklabel',['n=',num2str(sum(~isnan(x)))])
  return
end

a(1,1) = min(x(~g));
a(1,2) = max(x(~g));
a(1,3) = nanmean(x(~g));
a(1,4) = mode(x(~g));
a(1,4) = nanstd(x(~g));

a(2,1) = min(x(~~g));
a(2,2) = max(x(~~g));
a(2,3) = nanmean(x(~~g));
a(2,4) = nanstd(x(~~g));

cmap = [
  0.4660    0.6740    0.1880
  0.8500    0.3250    0.0980];
if exist('xname','var'), ylabel(xname); end
if numel(unique(x)) == 2
  hold on
  h1 = bar(0, mean(x(~g)),0.5, 'facecolor',cmap(1,:),'facealpha',0.5, 'linestyle','none');
  h2 = bar(1, mean(x(~~g)),0.5, 'facecolor',cmap(2,:),'facealpha',0.5, 'linestyle','none');
  hold off
  box on
  xlim([-0.5 1.5])
  [chi2,p,df] = chi2test(~~x,g+1);
  title(sprintf('{\\itx^2}(%d)=%.2f, {\\itp}=%.3f',df,chi2,p))
else
  cfg=struct('showdata',3,'groupcolor',cmap,'outliers',0,'maxwhisker',1.5);
  cat_plot_boxplot_sg({x(~g), x(~~g)}, cfg)
  [~,p,~,st] = ttest2(x(~g), x(~~g));
  title(sprintf('{\\itt}(%d)=%.2f, {\\itp}=%.3f',st.df,st.tstat,p))
end
xticklabel = {sprintf('Not %s (n=%d)',groupname,sum(~isnan(x(~g)))), ...
  sprintf('%s (n=%d)',groupname,sum(~isnan(x(~~g))))};
set(gca,'xticklabel',xticklabel, 'xticklabelrotation',90)
ylabel(xname)
end