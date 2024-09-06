function mytrf_plotweights(ax, nfeats, chanlocs, t, w)
zw = zscore(w,[],2);
for ifeat = 1:nfeats
  axespos(ax,8*(ifeat)+1)
  % size(w) = #feats x #lags x #chans x #folds x #subs
  plot(t, squeeze(mean(mean(w(ifeat,:,:,:,:),3),4)) )
  ylabel('\beta (all chan)'); xlabel('Time [ms]')
  title('Individual')
  xlim([t(1)+50 t(end)-50])
  
  axespos(ax,8*(ifeat)+2)
  errorplot(t, squeeze(mean(mean(w(ifeat,:,:,:),3),4))')
  toi = ((0)<t) & (t<(t(end)-50));
  [~,idx_t1] = max(squeeze(mean(mean(w(ifeat,toi,:,:),3),4)));
  [~,idx_t2] = min(squeeze(mean(mean(w(ifeat,toi,:,:),3),4)));
  idx_t1 = idx_t1 + find(toi,1,'first');
  idx_t2 = idx_t2 + find(toi,1,'first');
  hold on
  ylim0 = ylim;
  plot([t(idx_t1) t(idx_t1)], ylim0, 'r')
  plot([t(idx_t2) t(idx_t2)], ylim0, 'b')
  ylabel('Avg \beta (all chan) \pm SE'); xlabel('Time [ms]')
  title('Group')
  xlim([t(1)+50 t(end)-50])
  
  axespos(ax,8*(ifeat)+3)
  topoplot( squeeze(mean(mean(w(ifeat,idx_t1,:,:),4),5)), chanlocs, ...
    'colormap',flipud(brewermap(256,'RdBu')));
  colorbaro
  title(sprintf('mean \\beta: %.0f ms', t(idx_t1)), ...
    'backgroundcolor','r','color','w')
  
  axespos(ax,8*(ifeat)+4)
  topoplot( squeeze(mean(mean(w(ifeat,idx_t2,:,:),4),5)), chanlocs, ...
    'colormap',flipud(brewermap(256,'RdBu')));
  colorbaro
  title(sprintf('mean \\beta: %.0f ms', t(idx_t2)), ...
    'backgroundcolor','b','color','w')
  
  %% Z-scored?
  axespos(ax,8*(ifeat)+5)
  plot(t, squeeze(mean(mean(zw(ifeat,:,:,:,:),3),4)))
  ylabel('Z(\beta) (all chan)'); xlabel('Time [ms]')
  xlim([t(1)+50 t(end)-50])
  
  axespos(ax,8*(ifeat)+6)
  errorplot(t, squeeze(mean(mean(zw(ifeat,:,:,:,:),3),4))')
  [~,idx_t1] = max(squeeze(mean(mean(zw(ifeat,toi,:,:),3),4)));
  [~,idx_t2] = min(squeeze(mean(mean(zw(ifeat,toi,:,:),3),4)));
  idx_t1 = idx_t1 + find(toi,1,'first');
  idx_t2 = idx_t2 + find(toi,1,'first');
  hold on
  ylim0 = ylim;
  plot([t(idx_t1) t(idx_t1)], ylim0, 'm')
  plot([t(idx_t2) t(idx_t2)], ylim0, 'g')
  ylabel('Avg Z(\beta) (all chan) \pm SE'); xlabel('Time [ms]')
  xlim([t(1)+50 t(end)-50])
  
  axespos(ax,8*(ifeat)+7)
  topoplot( squeeze(mean(mean(zw(ifeat,idx_t1,:,:),4),5)), chanlocs, ...
    'colormap',flipud(brewermap(256,'PiYG')));
  colorbaro
  title(sprintf('mean Z(\\beta): %.0f ms', t(idx_t1)),...
    'backgroundcolor','m','color','k')
  
  axespos(ax,8*(ifeat)+8)
  topoplot( squeeze(mean(mean(zw(ifeat,idx_t2,:,:),4),5)), chanlocs, ...
    'colormap',flipud(brewermap(256,'PiYG')));
  colorbaro
  title(sprintf('mean Z(\\beta): %.0f ms', t(idx_t2)), ...
    'backgroundcolor','g','color','k')
end
end