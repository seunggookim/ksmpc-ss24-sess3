function stat = myft_ttest1_fdr(CFG, timelock)
% CFG may use:

fn_data=CFG.fn_data;
[p1,f1,e1]=myfileparts(fn_data);
prefix=[p1,'/',f1,'_-',num2str(CFG.prestim),'-',num2str(CFG.poststim),'s'];
if isfield(CFG,'desname'), prefix=[prefix,'_',CFG.desname]; end
if ~isfield(CFG,'correctm'), CFG.correctm='no'; end
prefix=[prefix,'_t1-',CFG.correctm];
correctm=[CFG.correctm,', alpha=',num2str(CFG.alpha)];
if ~isfield(CFG,'overwrite'), CFG.overwrite=0; end

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
    cfg.method='analytic';
    cfg.correctm='fdr';
    stat.(name) = ft_timelockstatistics(cfg, timelock.(name));
  end
  save(fn_stat,'stat')
else
  load(fn_stat,'stat')
end

%% visualize:
for k=1:numel(CFG.cond)
  name=CFG.cond(k).name;
  titlestr={CFG.fn_data, ['T-test1: ',name,', ',correctm]};
  myft_plot_stat(stat, name, titlestr)
  fn_png=[prefix,'_',name,'_STAT.png'];
  screen2png(fn_png,80,0)
  close(gcf)
end