function myft_checkIC1(comp,c)
%%
nTrl=numel(comp.trial);

ax=axeslayout([2 3]);
axespos(ax,1)
cfg=[];
cfg.component=c;
cfg.layout='easycapM24.mat';
cfg.comment='no';
cfg.markersymbol='.';
ft_topoplotIC(cfg, comp)
colormap(gca,flipud(brewermap(256,'spectral')))

if nTrl>1
  
  X=cat(3,comp.trial{:});
  axespos(ax,2)
  hold on
  mypcolor(comp.time{1},1:numel(comp.trial),squeeze(X(c,:,:))');
  caxis0=caxis;
  caxis([-max(abs(caxis0)) max(abs(caxis0))]./2);
  %colormap(gca,jet)
  colormap(gca,flipud(brewermap(256,'RdYlBu')))
  %colormap(gca,sgcolormap('BCWYR256'));
  ylabel('Trial index');
  xlabel('Time [s]');
  title('IC weights [AU]')
  
  axespos(ax,5)
  hold on
  errorplot(comp.time{1}, squeeze(X(c,:,:))')
  xlabel('Time [s]'); ylabel('IC wegith [AU]');
  grid on; box on;
  xlim([0 1])
  
  axespos(ax,3)
  scatter(squeeze(mean(X(c,:,:),2)),1:nTrl,'filled',...
    'markerfacecolor',[0 0 0],'markerFaceAlpha',0.2)
  xlabel('Mean weight [AU]'); ylabel('Trial index')
  grid on; box on;
  
  axespos(ax,6)
  scatter(squeeze(std(X(c,:,:),[],2)),1:nTrl,'filled',...
    'markerfacecolor',[0 0 0],'markerFaceAlpha',0.2)
  xlabel('Std weight [AU]'); ylabel('Trial index')
  grid on; box on;
  
  axespos(ax,4)
  [~,f]=pwelch(comp.time{1}(1,:),[],[],[],comp.fsample);
  PSD=zeros(nTrl,numel(f));
  for t=1:nTrl % 141 samples
    [PSD(t,:),f]=pwelch(comp.trial{t}(c,:),[],[],[],comp.fsample);
  end
  [pks,idx]=findpeaks(mean(PSD));
  [~,idx_idx]=max(pks);
  idx=idx(idx_idx);
  pk=pks(idx_idx);
  h=errorplot(f, PSD);
  xlabel('Frequency [Hz]'); ylabel('W-PSD [au^2/Hz]')
  set(gca,'xscale','lin','yscale','lin','XMinorGrid','on','Yminorgrid','off')
  hold on;
  scatter(f(idx),pk,'bx')
else
  axespos(ax,2)
  [PSD(1,:),f]=pwelch(comp.trial{1}(c,:),[],[],[],comp.fsample);
  [pks,idx]=findpeaks(PSD);
  [~,idx_idx]=max(pks);
  idx=idx(idx_idx);
  pk=pks(idx_idx);
  h=errorplot(f, PSD);
  xlabel('Frequency [Hz]'); ylabel('W-PSD [au^2/Hz]')
  set(gca,'xscale','lin','yscale','lin','XMinorGrid','on','Yminorgrid','off')
  hold on;
  scatter(f(idx),pk,'bx')
  
  ax2=axeslayout([2 1],[.1 .05 .1 .15],[.02, .3, .02, .02]);
  axespos(ax2,2)
  x=comp.trial{1}(c,:);
  plot(comp.time{1}, x)
  xlabel('Time [s]');ylabel('IC weights [AU]');
  title('Whole data')
  
  axespos(ax,6)
  [~,s]=max(x(401:end-400));
  s=s+400;
  plot(comp.time{1}, x)
  xlim([comp.time{1}(s-400) comp.time{1}(s+400)])
  xlabel('Time [s]');ylabel('IC weights [AU]');
  title('Zoomed around max')
  
  axespos(ax,3)
  [acor,lag]=xcorr(x,4*comp.fsample);
  plot(lag/comp.fsample,acor,'b')
  xlim([-4 4])
  xlabel('Lag [s]'); ylabel('Autocorr');
end

end