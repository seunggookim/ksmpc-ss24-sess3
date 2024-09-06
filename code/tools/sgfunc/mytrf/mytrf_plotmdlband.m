function h = mytrf_plotmdlband(mse_opt, mse_eva, r_opt, r_eva, ...
  lambdas, lambda_opt, chanlocs, t, w)

nfeats = size(w,1);
nsub = size(mse_opt,3);
nL = numel(lambdas);
l10L = log10(lambdas);
%%
ax = axeslayout([1+nfeats,8],[.3, .03, .15, .25]);
ax.x(3:8:end) = ax.x(3:8:end) - 0.03;
ax.x(4:8:end) = ax.x(4:8:end) - 0.04;
ax.x(7:8:end) = ax.x(7:8:end) - 0.03;
ax.x(8:8:end) = ax.x(8:8:end) - 0.04;
h = figure('position',[1 1 8*200 (1+nfeats)*200],'visible','off');

%% PLOT optimization & evaluation results:
axespos(ax,1)
% cmap = repmat(get(gca,'ColorOrder'),[ceil(nsub/7) 1]);
plot(squeeze(mean(mse_opt,1)))
ylabel('mse_{val} (all chan)'); xlabel('lambdaset')

axespos(ax,2); hold on
imagesc(l10L, l10L, reshape(squeeze(mean(mean(mse_opt,1),3)),[nL nL]))
colorbaro('southoutside')
axis square
scatter(log10(lambda_opt(:,1)),log10(lambda_opt(:,2)),'r.')
x = log10(geomean(lambda_opt));
scatter(x(1), x(2), 'go')
colormap(gca,brewermap(128,'Blues'))
xlabel('log10');ylabel('log10')
title('MSE_{val}')

axespos(ax,3)
x = squeeze(std(mean(mse_eva,1),[],3));
topoplot( x, chanlocs, 'colormap',brewermap(128,'Blues'));
caxis([min(x) max(x)]); colorbaro
title('std (mse_{test})')

axespos(ax,4)
x = squeeze(mean(mean(mse_eva,1),3));
topoplot( x, chanlocs, 'colormap',brewermap(128,'Blues'));
caxis([min(x) max(x)]); colorbaro
title('mean (mse_{test})')

axespos(ax,5)
plot(squeeze(mean(r_opt,1)))
ylabel('r_{val} (all chan)'); xlabel('lambdaset')

axespos(ax,6); hold on
imagesc(l10L, l10L, reshape(squeeze(mean(mean(r_opt,1),3)),[nL nL]))
axis square
x = log10(geomean(lambda_opt));
scatter(x(1), x(2), 'co')
x = log10(lambda_opt);
scatter(x(:,1), x(:,2), 'w.')
colormap(gca,brewermap(128,'Reds'))
title('r_{val}')

axespos(ax,7)
topoplot( squeeze(std(mean(r_eva,1),[],3)), chanlocs, ...
  'colormap',brewermap(128,'Reds'));
clim0 = caxis; caxis([0 clim0(2)]); colorbaro
title('std (r_{test})')

axespos(ax,8)
topoplot( squeeze(mean(mean(r_eva,1),3)), chanlocs, ...
  'colormap',brewermap(128,'Reds'));
clim0 = caxis; caxis([0 clim0(2)]); colorbaro
title('mean (r_{test})')

%% PLOT weights (this part is regardless of optimization)
mytrf_plotweights(ax, nfeats, chanlocs, t, w)


set(gcf,'color','w')
end