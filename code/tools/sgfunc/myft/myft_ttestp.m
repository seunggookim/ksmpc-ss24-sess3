function [stat,CFG] = myft_ttestp(CFG, timelock)
% CFG may use:

fn_data=CFG.fn_data;
[p1,f1,e1]=myfileparts(fn_data);
prefix=[p1,'/',f1,'_-',num2str(CFG.prestim),'-',num2str(CFG.poststim),'s'];
if isfield(CFG,'desname'), prefix=[prefix,'_',CFG.desname]; end
if ~isfield(CFG,'correctm'), CFG.correctm='no'; end
if ~isfield(CFG,'correcttail'), CFG.correcttail='prob'; end
switch CFG.correctm
  case {'fdr','no','Bonferroi'}
    CFG.method='analytic';
  case {'cluster'}
    CFG.method='montecarlo';
end
names=fieldnames(timelock);
if strcmp(CFG.correctm,'cluster')
  if ~isfield(CFG,'numrandomization'), CFG.numrandomization=2000; end
  if ~isfield(CFG,'clusteralpha'), CFG.clusteralpha=0.01; end
  if ~isfield(CFG,'neighbours')
    CFG.neighbours=ft_prepare_neighbours(struct('method','triangulation'), ...
      timelock.(names{1}));
  end
  prefix_correctm=['clus-cft',num2str(CFG.clusteralpha),...
    '-lat',num2str(CFG.latency(1)),'-',num2str(CFG.latency(2)), ...
    '-perm',num2str(CFG.numrandomization)];
else
  prefix_correctm=[CFG.correctm,...
    '-lat',num2str(CFG.latency(1)),'-',num2str(CFG.latency(2))];
end
prefix=[prefix,'_tp-',prefix_correctm];
if ~isfield(CFG,'overwrite'), CFG.overwrite=0; end
if ~isfield(CFG,'kont') && isfield(CFG,'kontrastmat')
  K=CFG.kontrastmat;
  CFG.kont=[];
  for k=1:size(K,1)
    CFG.kont(k).cond1=names{K(k,1)};
    CFG.kont(k).cond2=names{K(k,2)};
    CFG.kont(k).name=[CFG.kont(k).cond1,'_vs_',CFG.kont(k).cond2];
  end
end
fn_stat=[prefix,'_STAT.mat'];
if ~exist(fn_stat,'file') || CFG.overwrite
  stat=[];
  for k=1:numel(CFG.kont)
    cond1=CFG.kont(k).cond1;
    cond2=CFG.kont(k).cond2;
    kont=CFG.kont(k).name;
    
    cfg=CFG;
    n1=size(timelock.(cond1).trial,1);
    n2=size(timelock.(cond2).trial,1);
    if n1~=n2
      error('The trials need to be paired!!')
    end
    cfg.design=[ones(1,n1) 2*ones(1,n2)
      repmat(1:n1,[1 2])];
    cfg.statistic='depsamplesT';
    cfg.ivar=1;
    cfg.uvar=2;
    cfg.clusterthreshold='nonparametric_individual';
    stat.(kont) = ft_timelockstatistics(cfg, ...
      timelock.(cond1), timelock.(cond2));
    stat.(kont).df = n1-1;
  end
  save(fn_stat,'stat')
else
  load(fn_stat,'stat')
end

%% visualize:
for k=1:numel(CFG.kont)
  kont=CFG.kont(k).name;
  correctm=[prefix_correctm,', alpha=',num2str(CFG.alpha)];
  titlestr={CFG.fn_data, ['Paried T-test: ',kont,', ',correctm]};
  CFG = myft_plot_stat(stat, kont, titlestr, CFG);
  
  fn_png=[prefix,'_',kont,'_STAT.png'];
  screen2png(fn_png,80,0)
  close(gcf)
end