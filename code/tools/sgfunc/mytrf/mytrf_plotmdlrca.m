function h = mytrf_plotmdlrca(mse_opt, mse_eva, r_opt, r_eva, ...
  lambdas, lambda_opt, chanlocs, t, w, rca)

% size(w) = #feats x #lags x #chans x #folds x #subs
ncomps = size(w,3);
% nfeats = size(w,1);
nsub = size(mse_opt,3);
%%
ax = axeslayout([1+ncomps,4],[.3, .03, .15, .25]);
ax.x(3:8:end) = ax.x(3:8:end) - 0.03;
ax.x(4:8:end) = ax.x(4:8:end) - 0.04;
ax.x(7:8:end) = ax.x(7:8:end) - 0.03;
ax.x(8:8:end) = ax.x(8:8:end) - 0.04;
h = figure('position',[1 1 4*200 (1+ncomps)*200],'visible','off');

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
ylabel('mse_{val} (all comp)');xlabel('log10(\lambda)')
ylim([mine maxe])
title(sprintf('std log10(L)=%.2f', std(log10(lambda_opt))))

axespos(ax,2); hold on
err = squeeze(mean(mse_opt,1));
errorplot(log10(lambdas), err')
ylim0 = ylim;
plot([1 1]*mean(log10(lambda_opt)), ylim0, 'color','r')
ylim(ylim0)
ylabel('mse_{val} (all comp)');xlabel('log10(\lambda)')
title(sprintf('mean log10(L)=%.2f', mean(log10(lambda_opt))))

axespos(ax,3)
hold on
for isub = 1:nsub
  r = squeeze(mean(r_opt(:,:,isub),1));
  plot(log10(lambdas), r, 'color',cmap(isub,:))
  plot([1 1]*log10(lambda_opt(isub)), [min(r) max(r)],'color',cmap(isub,:))
end
ylabel('r_{val} (all chan)'); xlabel('log10(\lambda)')

axespos(ax,4); hold on
r = squeeze(mean(r_opt,1));
errorplot(log10(lambdas), r')
ylabel('r_{val} (all comp)'); xlabel('log10(\lambda)')
plot([1 1]*mean(log10(lambda_opt)), [min(mean(r,2)) max(mean(r,2))], ...
  'color','r')

%% PLOT weights (this part is regardless of optimization)
mytrf_plotweightsrca(ax, ncomps, chanlocs, t, w, rca, r_eva)

set(gcf,'color','w')
end
