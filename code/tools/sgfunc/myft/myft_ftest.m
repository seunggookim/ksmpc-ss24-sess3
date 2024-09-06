function [stat,CFG] = myft_ftest(CFG, timelock)
% CFG may use:

fn_data=CFG.fn_data;
[p1,f1,e1]=myfileparts(fn_data);
prefix=[p1,'/',f1,'_-',num2str(CFG.prestim),'-',num2str(CFG.poststim),'s'];
if isfield(CFG,'desname'), prefix=[prefix,'_',CFG.desname]; end
if ~isfield(CFG,'correctm'), CFG.correctm='no'; end
%if ~isfield(CFG,'correcttail'), CFG.correcttail='prob'; end
CFG.tail=1;
CFG.correcttail='no';
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
prefix=[prefix,'_f-',prefix_correctm];
if ~isfield(CFG,'overwrite'), CFG.overwrite=0; end
%%
CFG.kont=[];
CFG.kont(1).name=CFG.desname;

fn_stat=[prefix,'_STAT.mat'];
if ~exist(fn_stat,'file') || CFG.overwrite
  stat=[];
  for k=1:numel(CFG.kont)
    kont=CFG.kont(k).name;
    cfg=CFG;
    cfg.design=[];
    for j=1:numel(names)
      cfg.design=[cfg.design j*ones(1,size(timelock.(names{j}).(CFG.parameter),1))];
    end
    cfg.statistic='indepsamplesF';
    cfg.ivar=1;
    cmd=['stat.(kont) = ft_timelockstatistics(cfg, '];
    for j=1:numel(names)
      cmd=[cmd, 'timelock.',names{j},','];
    end
    cmd=[cmd(1:end-1),');'];
    eval(cmd)
  end
  save(fn_stat,'stat')
else
  load(fn_stat,'stat')
end

%% visualize:
for k=1:numel(CFG.kont)
  kont=CFG.kont(k).name;
  correctm=[prefix_correctm,', alpha=',num2str(CFG.alpha)];
  titlestr={CFG.fn_data, ['F-test: ',kont,', ',correctm]};
  CFG = myft_plot_stat(stat, kont, titlestr, CFG);
  
  fn_png=[prefix,'_',kont,'_STAT.png'];
  screen2png(fn_png,80,0)
  close(gcf)
end