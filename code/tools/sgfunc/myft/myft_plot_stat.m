function CFG=myft_plot_stat(stat, name, titlestr, CFG)

if strcmp(stat.change.cfg.statistic(end),'F')
  isF=1; statname='F';
else
  isF=0; statname='T';
end

figure('position',[  -759   505   740   450])
ax1=axeslayout([4 1],[.1 0 .2 .15],[.02, .02, .02, .01]);
axespos(ax1,1)
maskedT=stat.(name).stat.*stat.(name).mask;
mypcolor(stat.(name).time, 1:157, maskedT)
cmap=flipud(brewermap(256,'RdBu'));
if isF
  cmap=flipud(brewermap(256*2,'RdBu'));
  cmap=cmap(257:end,:);
end
colormap(gca, cmap);
c0=caxis;
c0=[-max(abs(c0)) max(abs(c0))];
if isF, c0=[0 c0(2)]; end
caxis(c0)
title(titlestr,'interpreter','none', 'fontsize',8);
hcb=colorbar;
if isF
  if isfield(stat.(name),'dfnum') && isfield(stat.(name),'dfdenom')
    title(hcb,[statname,'(',num2str(stat.(name).dfnum),',',...
      num2str(stat.(name).dfdenom),')'])
  else
    title(hcb,[statname,'(nan)'])
  end
else
  if isfield(stat.(name),'df'), df=stat.(name).df(1); else df=nan; end
  title(hcb,[statname,'(',num2str(df),')'])
end
ylabel('Channel');

ax3=axeslayout([4 1],[.1 .099 .1 .15],[.02, .02, .02, .01]);
axespos(ax3,2)
plot(stat.(name).time, stat.(name).stat, 'color',[0 0 0 .1])
set(gca,'xlim',[stat.(name).time(1) stat.(name).time(end)])
xlabel('Time [s]'); ylabel([statname,'-stat'])
y0=ylim;
y0=[-max(abs(ylim)) max(abs(ylim))];
if isF, y0=[0 y0(2)]; end

if isfield(stat.(name),'posclusters') || isfield(stat.(name),'negclusters')
  % find time interval of significant clusters:
  XLIM={};PEAK=[];
  % find significant positive clusters
  if isfield(stat.(name),'posclusters')
    clus=stat.(name).posclusters;
    idx_clus=find([clus.prob]<stat.(name).cfg.alpha);
    for j=1:numel(idx_clus)
      mask=stat.(name).posclusterslabelmat==idx_clus(j);
      trms=rms(denan(stat.(name).stat.*mask),1);
      [~,idx_peak]=max(trms);
      PEAK=[PEAK stat.(name).time(idx_peak)];
      ind_time=any(mask,1);
      xlim=stat.(name).time(ind_time);
      xlim=xlim([1 end]);
      XLIM=[XLIM xlim];
    end
  end
  % find significant negative clusters
  if isfield(stat.(name),'negclusters')
    clus=stat.(name).negclusters;
    idx_clus=find([clus.prob]<stat.(name).cfg.alpha);
    for j=1:numel(idx_clus)
      %ind_time=any(stat.(name).negclusterslabelmat==idx_clus(j),1);
      mask=stat.(name).negclusterslabelmat==idx_clus(j);
      trms=rms(denan(stat.(name).stat.*mask),1);
      [~,idx_peak]=max(trms);
      PEAK=[PEAK stat.(name).time(idx_peak)];
      ind_time=any(mask,1);
      xlim=stat.(name).time(ind_time);
      xlim=xlim([1 end]);
      XLIM=[XLIM xlim];
    end
  end
  nClus=numel(PEAK);
  PEAK=sort(PEAK); % reorder time points
else
  % find RMS peaks
  y=rms(denan(maskedT));
  [pks,idx_pks]=findpeaks(y);
  PEAK=stat.(name).time(idx_pks);
  nClus=numel(PEAK);
end

% mark timepoint on the plot
linecolor=brewermap(nClus,'Dark2');
for j=1:nClus
  hold on
  line([PEAK(j) PEAK(j)]',y0, 'color',linecolor(j,:))
end

% draw topolot
dt=stat.(name).time(2)-stat.(name).time(1);
ax2=axeslayout([2 nClus+1],[0 0 0.2 0.2],[.02, .02, .02, .01]);
for j=1:nClus
  axespos(ax2,j+(nClus+1))
  toi=[PEAK(j)-dt PEAK(j)+dt];
  cfg=struct('xlim',toi,'highlightcolor','k','zlim',c0);
  myft_topo(cfg, stat.(name))
  colormap(gca, flipud(brewermap(256,'RdBu')));
  title(sprintf('%.3f s',PEAK(j)),'fontsize',8, 'color',linecolor(j,:));
end
axespos(ax2,2*(nClus+1))
cfg=struct('xlim',CFG.latency,'style','blank','highlightcolor','k');
cfg=myft_topo(cfg, stat.(name));
if ~isfield(cfg,'highlightchannel'), cfg.highlightchannel=[]; end
title(sprintf('%.3f-%.3f\nN=%i',CFG.latency,...
  numel(cfg.highlightchannel)),...
  'fontsize',8, 'color','k');
end