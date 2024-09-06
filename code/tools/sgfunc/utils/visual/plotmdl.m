function h = plotmdl(mdl,v)
if exist('v','var')
  h = plotAdded(mdl,v);
  jdx = find(ismember(mdl.CoefficientNames, v));
  if numel(mdl.VariableNames)>2
    prefix = 'Adj. ';
  else
    prefix = '';
  end
else
  h = plot(mdl);
  prefix = '';
  idx = 1;
end

if numel(mdl.CoefficientNames) > 1
  title(sprintf('EDF=%i, Adj. R^2=%.3f', mdl.DFE, mdl.Rsquared.Adjusted))
  xlabel(sprintf('%s%s: est=%.3f, p=%.0i',...
    prefix, v, mdl.Coefficients.Estimate(jdx), ...
    mdl.Coefficients.pValue(jdx)), 'interp','none')
  legend off
  x = h(1).XData;
  y = h(1).YData;
  h(1).Visible = 'off';
  hold on
  n = mdl.NumObservations;
  alpha = rescale(-log(n),0,0.5,'inputmin',-19,'inputmax',0);
  h(1) = scatter(x,y,20, 'MarkerFaceColor',[0 .6 .6],'MarkerFaceAlpha',alpha,...
    'MarkerEdgeColor','none');
  h(2).Color = [1 0 .2];
  if exist('v','var')
    x = h(3).XData;
    y = h(3).YData;
    h(3).Visible = 'off';
  else
    x = [h(3).XData, nan, h(4).XData];
    y = [h(3).YData, nan, h(4).YData];
    h(3).Visible = 'off';
    h(4).Visible = 'off';
  end
  idx = find(isnan(x));
  h(5) = patch([x(1:idx-1) fliplr(x(idx+1:end))],...
    [y(1:idx-1) fliplr(y(idx+1:end))], [1 0 0.2], 'linestyle','none', ...
    'facealpha',0.3);
elseif numel(mdl.CoefficientNames) == 1
  % point-estimate: what R^2?
  histogram(mdl.Residuals.Raw + mdl.Coefficients.Estimate)
  xline(mdl.Coefficients.Estimate,'r')
  title(sprintf('T(%i)=%.2f, RMSE=%.2f', ...
    mdl.DFE, mdl.Coefficients.tStat, mdl.RMSE))
  xlabel(sprintf('%s%s: est=%.3f, p=%.0i',...
    prefix, mdl.VariableNames{idx}, mdl.Coefficients.Estimate, ...
    mdl.Coefficients.pValue),'interp','none')
  ylabel('#samples')

else
  error('what else?')
end

hold off
end
