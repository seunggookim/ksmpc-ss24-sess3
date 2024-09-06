function [epoched,CFG] = myft_epoch(CFG)
% CFG requires:
%  .fn_data  '1xN' filename for preprocessed data (.mat)
%  .prestim  [1x1] time interval before stim onset [sec]
%  .poststim [1x1] time interval after stim onset [sec]
%  .cond     (1xK) condition structure
%  .cond(k).name '1xN' condition name
%  .cond(k).keyval {1xP} key-value pairs

%%
fn_data = CFG.fn_data;
[p1,f1,e1] = myfileparts(fn_data);

% define trials
ls(fn_data)
load(fn_data,'data'); % read data

% select MEG channels only
data = ft_selectdata(struct('channel','AG*'), data);

if ~isfield(CFG,'fn_tab')
  fn_tab=[p1,'/onset_stim.xls'];
else
  fn_tab=[p1,'/',CFG.fn_tab];
end
ls(fn_tab)
tab = readtable(fn_tab); % read stim onsets
t = tab.onset_sec; 
trl = [t-CFG.prestim t+CFG.poststim t*0-CFG.prestim]; % all trials
trl = round(trl*data.fsample);
idx_tab = ~(trl(:,1)<1 | trl(:,2)>numel(data.time{1}));
trl = trl(idx_tab,:);

% epoch ALL events
cfg = struct('trl',trl);%, 'minlength',CFG.prestim+CFG.poststim);
ft_warning off; ft_debug off
data = ft_redefinetrial(cfg, data);
% warning('... because of resampling after loading the raw data')

% find trials to reject
cfg = CFG;
cfg.trl = trl;
cfg.zthres_dvnt = 10; % otherwise too many trials rejected
[isJmpyTrl,isDvnTrl]  = myft_reject_trials(cfg, data);

% define trials by conditions
epoched = [];
for k = 1:numel(CFG.cond)
  nKeyvalpair = numel(CFG.cond(k).keyval)/2;
  idx = [~isJmpyTrl ~isDvnTrl]; % rejecting "bad" trials
  for p = 1:nKeyvalpair
    trgkey = CFG.cond(k).keyval{1+(p-1)*2};
    trgval = CFG.cond(k).keyval{p*2};
    trgkey = tab.(trgkey);
    idx = [idx trgkey(idx_tab)==trgval];
  end
  CFG.cond(k).idx_trial = find(all(idx,2)); % is it always "all"?
  nTrl = numel(CFG.cond(k).idx_trial);
  fprintf('found %d trials for "%s"\n', nTrl, CFG.cond(k).name)
  
  % select trials by conditions
  cfg=struct('trials',CFG.cond(k).idx_trial);
  epoched.(CFG.cond(k).name) = ft_selectdata(cfg, data);
end

end