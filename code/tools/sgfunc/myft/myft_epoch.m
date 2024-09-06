function [epoched,CFG] = myft_epoch(CFG)
% CFG requires:
%  .fn_data  '1xN' filename for preprocessed data (.mat)
%  .prestim  [1x1] time interval before stim onset [sec]
%  .poststim [1x1] time interval after stim onset [sec]
%  .cond     (1xK) condition structure
%  .cond(k).name '1xN' condition name
%  .cond(k).keyval {1xP} key-value pairs

% load data
fn_data = CFG.fn_data;
ls(fn_data)
load(fn_data,'data'); % read data
% data = ft_selectdata(struct('channel','AG*'), data); % select MEG channels only

% filtering
if sum(contains(fieldnames(CFG),'filter'))
  data = ft_preprocessing(CFG, data);
end

% find stimuli of interest
[p1,~,~] = myfileparts(fn_data);
if ~isfield(CFG,'fn_tab')
  fn_tab=[p1,'/onset_stim.xls'];
else
  fn_tab=[p1,'/',CFG.fn_tab];
end
ls(fn_tab)
tab = readtable(fn_tab); % read stim onsets
trl_cond = [];
for k = 1:numel(CFG.cond) % find stimuli that satisfy key-value conditions
  nKeyvalpair = numel(CFG.cond(k).keyval)/2;
  idx = [];
  for p = 1:nKeyvalpair
    trgkey = CFG.cond(k).keyval{1+(p-1)*2};
    trgval = CFG.cond(k).keyval{p*2};
    trgkey = tab.(trgkey);
    idx = [idx trgkey==trgval];
  end
  CFG.cond(k).idx_tab = find(all(idx,2));
  trl_cond=[trl_cond; ones(size(CFG.cond(k).idx_tab))*k];
end
t = tab.onset_sec(cat(1,CFG.cond.idx_tab)); % onset time of all stim of interest
trl = [t-CFG.prestim t+CFG.poststim t*0-CFG.prestim]; % all trials of interest [sec]
trl = round(trl*data.fsample); % [sample]
trl = [trl trl_cond]; % [start_smp end_smp offset_smp condIndex]
ind_tab = ~(trl(:,1)<1 | trl(:,2)>numel(data.time{1})); % exclude trials with sample indices outside of the data
trl = trl(ind_tab,:);

% epoch ALL events
cfg = struct('trl',trl);%, 'minlength',CFG.prestim+CFG.poststim);
ft_warning off; ft_debug off
data = ft_redefinetrial(cfg, data);
% warning('... because of resampling after loading the raw data')

% find trials to reject
cfg = CFG;
cfg.trl = trl;
cfg.zthres_dvnt = 10; % otherwise too many trials rejected
[isJmpyTrl,isDvnTrl] = myft_reject_trials(cfg, data);

% define trials by conditions
epoched = [];
for k = 1:numel(CFG.cond)
  CFG.cond(k).idx_trl = find(trl(:,4) == k & ~isJmpyTrl & ~isDvnTrl); % rejecting "bad" trials
  nTrl = numel(CFG.cond(k).idx_trl);
  fprintf('found %d trials for "%s"\n', nTrl, CFG.cond(k).name)
  
  % select trials by conditions
  cfg=struct('trials',CFG.cond(k).idx_trl);
  epoched.(CFG.cond(k).name) = ft_selectdata(cfg, data);
end

end