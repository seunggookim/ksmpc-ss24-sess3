function stat = myft_ttest2(CFG, timelock)
% CFG may use:

fn_data=CFG.fn_data;
[p1,f1,e1]=myfileparts(fn_data);
prefix=[p1,'/',f1,'_-',num2str(CFG.prestim),'-',num2str(CFG.poststim),'s'];
if isfield(CFG,'desname'), prefix=[prefix,'_',CFG.desname]; end
if ~isfield(CFG,'correctm'), CFG.correctm='fdr'; end
prefix=[prefix,'_t2-',CFG.correctm];
correctm=[CFG.correctm,', alpha=',num2str(CFG.alpha)];
if ~isfield(CFG,'overwrite'), CFG.overwrite=0; end

fn_stat=[prefix,'_STAT.mat'];
if ~exist(fn_stat,'file') || CFG.overwrite
  stat=[];
  for k=1:numel(CFG.kont)
    cond1=CFG.kont(k).cond1;
    cond2=CFG.kont(k).cond2;
    kont=[cond1,'_vs_',cond2];
    
    cfg=CFG;
    cfg.design=[ones(1,size(timelock.(cond1).trial,1)) ...
      2*ones(1,size(timelock.(cond2).trial,1))];
    cfg.statistic='indepsamplesT';
    cfg.ivar=1;
    cfg.method='analytic';
    stat.(kont) = ft_timelockstatistics(cfg, timelock.(cond1), timelock.(cond2));
  end
  save(fn_stat,'stat')
else
  load(fn_stat,'stat')
end

%% visualize:
for k=1:numel(CFG.cond)
  kont=CFG.cond(k).name;
  titlestr={CFG.fn_data, ['T-test2: ',kont,', ',correctm]};
  myft_plot_stat(stat, kont, titlestr)
  fn_png=[prefix,'_',kont,'_STAT.png'];
  screen2png(fn_png,80,0)
  close(gcf)
end