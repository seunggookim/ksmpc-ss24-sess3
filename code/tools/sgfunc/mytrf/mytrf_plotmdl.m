function h = mytrf_plotmdl(mse_opt, mse_eva, r_opt, r_eva, ...
  lambdas, lambda_opt, chanlocs, t, w)

nfeats = size(w,1);
nsub = size(mse_opt,3);
%%
ax = axeslayout([1+nfeats,8],[.3, .03, .15, .25]);
ax.x(3:8:end) = ax.x(3:8:end) - 0.03;
ax.x(4:8:end) = ax.x(4:8:end) - 0.04;
ax.x(7:8:end) = ax.x(7:8:end) - 0.03;
ax.x(8:8:end) = ax.x(8:8:end) - 0.04;
h = figure('position',[1 1 8*200 (1+nfeats)*200],'visible','on');

%% PLOT optimization & evaluation results:
axespos(ax,1)
cmap = repmat(get(gca,'ColorOrder'),[ceil(nsub/7) 1]);
hold on
mine = 1; maxe = 1;
for isub = 1:nsub
  err = squeeze(mean(mse_opt(:,:,isub),1));
  plot(log10(lambdas), err, 'color',cmap(isub,:))
  plot([1 1]*log10(lambda_opt(isub)), [0 min(err)],'color',cmap(isub,:))
  mine = min(mine, min(err));
  maxe = max(maxe, max(err)); 
end
ylabel('mse_{val} (all chan)');xlabel('log10(\lambda)')
ylim([mine maxe])
title(sprintf('std log10(L)=%.2f', std(log10(lambda_opt))))

axespos(ax,2); hold on
err = squeeze(mean(mse_opt,1));
errorplot(log10(lambdas), err')
ylim0 = ylim;
plot([1 1]*mean(log10(lambda_opt)), ylim0, 'color','r')
ylim(ylim0)
ylabel('mse_{val} (all chan)');xlabel('log10(\lambda)')
title(sprintf('mean log10(L)=%.2f', mean(log10(lambda_opt))))

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
hold on
for isub = 1:nsub
  r = squeeze(mean(r_opt(:,:,isub),1));
  plot(log10(lambdas), r, 'color',cmap(isub,:))
  plot([1 1]*log10(lambda_opt(isub)), [min(r) max(r)],'color',cmap(isub,:))
end
ylabel('r_{val} (all chan)'); xlabel('log10(\lambda)')

axespos(ax,6); hold on
r = squeeze(mean(r_opt,1));
errorplot(log10(lambdas), r')
ylabel('r_{val} (all chan)'); xlabel('log10(\lambda)')
plot([1 1]*mean(log10(lambda_opt)), [min(mean(r,2)) max(mean(r,2))], ...
  'color','r')

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