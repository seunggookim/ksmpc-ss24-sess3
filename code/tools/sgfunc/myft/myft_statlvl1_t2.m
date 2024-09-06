function stat = myft_statlvl1_t2(CFG, timelock)
% CFG may use:
% (.clusteralpha=0.01)
% (.numrandomization=2000)

fn_data=CFG.fn_data;
[p1,f1,e1]=myfileparts(fn_data);
prefix=[p1,'/',f1,'_-',num2str(CFG.prestim),'-',num2str(CFG.poststim),'s'];
prefix=[prefix,'_t2-clus'];
if ~isfield(CFG,'overwrite'), CFG.overwrite=0; end
if ~isfield(CFG,'numrandomization'), CFG.numrandomization=2000; end
if ~isfield(CFG,'clusteralpha'), CFG.clusteralpha=0.01; end

fn_stat=[prefix,'_STAT.mat'];
if ~exist(fn_stat,'file') || CFG.overwrite
  stat=[];
  for k=1:numel(CFG.cntrst)
    cntrstname=CFG.cntrst(k).name;
    
    cfg=CFG;
    cfg.design=ones(1,size(timelock.(name).trial,1));
    cfg.statistic='indepsamplesT';
    cfg.ivar=1; % row number for ivar in the design matrix
    cfg.method='montecarlo';
    cfg.correctm='cluster';
    %cfg.clusteralpha=CFG.clusteralpha;
    cfg.neighbours=ft_prepare_neighbours(struct('method','triangulation'), ...
      timelock.(name));
    %cfg.numrandomization=CFG.numrandomization;
    stat.(cntrstname) = ft_timelockstatistics(cfg, timelock.(name1), timelock.(name2));
  end
  save(fn_stat,'stat')
else
  load(fn_stat,'stat')
end

%% visualize:
for k=1:numel(CFG.cond)
  name=CFG.cond(k).name;
  
  figure('position',[  -759   505   740   450])
  ax1=axeslayout([2 1],[.1 0 .15 .15]);
  axespos(ax1,1)
  mypcolor(stat.(name).time, 1:157, stat.(name).stat.*stat.(name).mask)
  colormap(gca, flipud(brewermap(256,'RdBu')));
  c0=caxis;
  c0=[-max(abs(c0)) max(abs(c0))];
  caxis(c0)
  df=median(timelock.(name).dof(:));
  title({CFG.fn_data, ['T-test1: ',name,...
    ', CBPT',...
    '(cft=',num2str(stat.(name).cfg.clusteralpha),', '...
    '#rand=',num2str(stat.(name).cfg.numrandomization),')']},...
    'interpreter','none', 'fontsize',8);
  hcb=colorbar;
  title(hcb,['T(',num2str(df),')'])
  xlabel('Time [s]'); ylabel('Channel');
  
  % find time interval of significant clusters:
  XLIM={};PEAK=[];
  % find significant positive clusters
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
  % find significant negative clusters
  clus=stat.(name).negclusters;
  idx_clus=find([clus.prob]<stat.(name).cfg.alpha);
  for j=1:numel(idx_clus)
    %ind_time=any(stat.(name).negclusterslabelmat==idx_clus(j),1);
    mask=stat.(name).negclusterslabelmat==idx_clus(j);
    trms=rms(denan(stat.(name).stat.*mask),1);
    [~,idx_peak]=max(trms);
    PEAK=[PEAK stat.(name).time(idx_peak)];
    
    xlim=stat.(name).time(ind_time);
    xlim=xlim([1 end]);
    XLIM=[XLIM xlim];
  end
  nClus=numel(PEAK);
  PEAK=sort(PEAK); % reorder time points
  
  % mark timepoint on the plot
  linecolor=brewermap(nClus,'Dark2');
  for j=1:nClus
    hold on
    line([PEAK(j) PEAK(j)]',[1 157]', 'color',linecolor(j,:))
  end
  
  % draw topolot
  dt=stat.alltones.time(2)-stat.alltones.time(1);
  ax2=axeslayout([2 nClus],[0 0 0.1 0.1]);
  for j=1:nClus
    axespos(ax2,j+nClus)
    toi=[PEAK(j)-dt PEAK(j)+dt];
    cfg=struct('xlim',toi,'highlightcolor','k','zlim',c0);
    myft_topo(cfg, stat.(name))
    colormap(gca, flipud(brewermap(256,'RdBu')));
    title([num2str(PEAK(j)),' s'],'fontsize',8, 'color',linecolor(j,:));
  end
  fn_png=[prefix,'_',name,'_STAT.png'];
  screen2png(fn_png,150)
  close(gcf)
end