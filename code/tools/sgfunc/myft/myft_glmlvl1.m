function stat = myft_glmlvl1(CFG, timelock)
% CFG may use:
% (.clusteralpha=0.01)
% (.numrandomization=2000)

fn_data=CFG.fn_data;
[p1,f1,e1]=myfileparts(fn_data);
prefix=[p1,'/',f1,'_-',num2str(CFG.prestim),'-',num2str(CFG.poststim),'s'];
if ~isfield(CFG,'correctm'), CFG.correctm='cluster'; end
if ~isfield(CFG,'correcttail'), CFG.correcttail='prob'; end
if strcmp(CFG.correctm,'cluster')
  prefix=[prefix,'_t1-clus-cft',num2str(CFG.clusteralpha),...
  '-lat',num2str(CFG.latency(1)),'-',num2str(CFG.latency(2)), ...
  '-perm',num2str(CFG.numrandomization)];
else
  prefix=[prefix,'_t1-',CFG.correctm];
end
if ~isfield(CFG,'overwrite'), CFG.overwrite=0; end
if ~isfield(CFG,'numrandomization'), CFG.numrandomization=2000; end
if ~isfield(CFG,'clusteralpha'), CFG.clusteralpha=0.01; end

fn_stat=[prefix,'_STAT.mat'];
if ~exist(fn_stat,'file') || CFG.overwrite
  stat=[];
  for k=1:numel(CFG.cond)
    name=CFG.cond(k).name;
    cfg=CFG;
    cfg.design=ones(1,size(timelock.(name).trial,1));
    cfg.statistic='onesamplesT';
    cfg.resampling='swapsign';
    cfg.ivar=1;
    cfg.method='montecarlo';
    cfg.correctm='cluster';
    cfg.neighbours=ft_prepare_neighbours(struct('method','triangulation'), ...
      timelock.(name));
    stat.(name) = ft_timelockstatistics(cfg, timelock.(name));
  end
  save(fn_stat,'stat')
else
  load(fn_stat,'stat')
end

%% visualize:
for k=1:numel(CFG.cond)
  name=CFG.cond(k).name;
  titlestr={CFG.fn_data, ['T-test1: ',name,...
    ', CBPT',...
    '(cft=',num2str(stat.(name).cfg.clusteralpha),', '...
    '#rand=',num2str(stat.(name).cfg.numrandomization),')']};
  myft_plot_stat(stat, name, titlestr, CFG)
  fn_png=[prefix,'_',name,'_STAT.png'];
  screen2png(fn_png,80)
  close(gcf)
end