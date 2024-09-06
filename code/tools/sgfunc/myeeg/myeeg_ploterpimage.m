function fn_png = myeeg_ploterpimage(eeg, cfg)

eeg = pop_eegfiltnew(eeg,'locutoff',0.5, 'hicutoff',20);
if ~exist('cfg','var')
  cfg = [];
end
if ~isfield(cfg,'events')
  cfg.events = unique({eeg.event(contains({eeg.event.code},'Stimulus')).type});
end
eeg = pop_epoch(eeg, cfg.events, [-0.5 6.5]);
if ~isfield(cfg,'ichan')
  cfg.ichan = 3;
end
if ~isfield(cfg,'titleprefix')
  cfg.titleprefix = '';
end
figure('visible','off')
pop_erpimage(eeg, 1, cfg.ichan, [], ...
  sprintf('%s%s, [0.5,20] Hz, (%i/%i) trials', ...
  cfg.titleprefix, eeg.chanlocs(cfg.ichan).labels, eeg.trials, ...
  eeg.p.n_trials_total),...
  0,1,{},[],'' ,'yerplabel','\muV','erp','on','cbar','on')
[~,f1,~] = fileparts(eeg.filename);
colormap(flipud(brewermap(256,'RdBu')));
h = get(gcf,'children');
set(h(2),'color','w')
hh = get(h(4),'children');
hh.CDataMapping = 'scaled';
set(gcf,'color','w')
fn_png = fullfile(eeg.filepath, [f1,'_erpimage.png']);
export_fig(fn_png,'-r120')
close(gcf)
end
