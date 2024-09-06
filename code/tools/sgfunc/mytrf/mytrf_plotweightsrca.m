function mytrf_plotweightsrca(ax, ncomps, chanlocs, t, w, rca, r_eva)
zw = zscore(w,[],2);
for icomp = 1:ncomps
  axespos(ax,4*(icomp)+1)
  % size(w) = #feats x #lags x #chans x #folds x #subs
  hold on
  cmap = get(gca,'colororder');
  for ifeat = 1:size(w,1)
    plot(t, squeeze(mean(w(ifeat,:,icomp,:,:),4)),'color',[cmap(ifeat,:) 0.5])
  end
  ylabel('\beta (all feats)'); xlabel('Time [ms]')
  title('Individual')
  xlim([t(1)+50 t(end)-50])
  
  axespos(ax,4*(icomp)+2)
  hold on
  hs = []; featindex = {};
  for ifeat = 1:size(w,1)
    h = errorplot(t, squeeze(mean(w(ifeat,:,icomp,:,:),4))');
    h.patch.FaceColor = cmap(ifeat,:);
    h.patch.FaceAlpha = 0.5;
    h.mainLine.Color = cmap(ifeat,:);
    h.edge(1).Color = cmap(ifeat,:);
    h.edge(2).Color = cmap(ifeat,:);
    featindex = [featindex ['Feat',num2str(ifeat)]];
    hs = [hs h.mainLine];
  end
  legend(hs, featindex)
  ylabel('Avg \beta (all chan) \pm SE'); xlabel('Time [ms]')
  title('Group')
  xlim([t(1)+50 t(end)-50])
  
  axespos(ax,4*(icomp)+3)
  topoplot( rca.A(:,icomp), chanlocs, ...
    'colormap',flipud(brewermap(256,'RdBu')))
%   colorbaro
  title('Fwd Mtx [A]')
  
  axespos(ax,4*(icomp)+4)
  % size(r_eva) = #chans x #folds x #subs
  meanr_eva = squeeze(mean(r_eva(icomp,:,:),2));
  hold on
  hs = scatter(meanr_eva*0+1+(rand(size(meanr_eva))-0.5)*0.2, meanr_eva,'.');
  hb = boxplot(meanr_eva,'colors',cmap);
  hs.CData = cmap(4,:);
  set(hb(contains(get(hb,'Tag'),'Outliers')), 'MarkerEdgeColor',cmap(3,:));
  ylabel('r')
end