function ica = myeeg_icaafterclean_sesICA(fnames_run, cfg)
% ica = myeeg_icaafterclean_sesICA(fnames_run, cfg)

if ~exist('cfg','var')
  cfg =[];
  cfg.HP = 1;         % cut-off frequency high-pass filter [Hz] only for ICA
  cfg.LP = 40;        % cut-off frequency low-pass filter [Hz] only for ICA
  cfg.SRATE = 250;    % downsample data for ICA
  cfg.HP_ord = [];   % high-pass filter order depends on sampling rate
  cfg.LP_ord = [];   % low-pass filter order depends on sampling rate
  cfg.PRUNE = 3;      % artifact rejection threshold in SD (for ICA only)
  cfg.PCA = 1;        % choose PCA option for ICA
  cfg.PCADIMS = 50;   % PCA dimension if PCA option is true
end

[p1,~,~] = fileparts(fnames_run{1});
fn_ica = fullfile(p1,'ica.mat');
if isfile(fn_ica)
  load(fn_ica, 'ica')
  return
end
fn_log = fullfile(p1,'ica.log');
if isfile(fn_log)
  delete(fn_log)
end
diary(fn_log)

% load rawdata:
for ifile = 1:numel(fnames_run)
  fn_run = fnames_run{ifile};
  eeg1 = pop_loadset(fn_run);
  
  chanlocs = eeg1.chanlocs;
  % Aggressively clean the raw data (only for ICA)
  eeg1 = pop_clean_rawdata(eeg1, ...
    'FlatlineCriterion',5, 'ChannelCriterion',0.8, 'LineNoiseCriterion',4,...
    'Highpass',[.25 .75], 'BurstCriterion',20, 'WindowCriterion',.25, ...
    'BurstRejection','on', 'Distance','Euclidian',...
    'WindowCriterionTolerances',[-Inf 7]);
  eeg1 = pop_interp(eeg1, chanlocs, 'spherical');
  
  % apply low pass filter (only for ICA)
  eeg1 = pop_eegfiltnew(eeg1, 'hicutoff',cfg.LP);
  % downsample data - may be an optional step
  eeg1 = pop_resample(eeg1, cfg.SRATE);
  % apply high pass filter (only for ICA)
  eeg1 = pop_eegfiltnew(eeg1, 'locutoff',cfg.HP);
  
  % create dummy events and epoch data to these dummy events
  eeg1 = eeg_regepochs(eeg1, 'recurrence', 1, 'eventtype', '999');
  eeg1 = eeg_checkset(eeg1, 'eventconsistency');
  
  % remove epochs with artefacts to improve ICA training
  eeg1 = pop_jointprob(eeg1, 1, [1:size(eeg1.data,1)], ...
    cfg.PRUNE, cfg.PRUNE, 0, 1, 0);
  
  if ifile == 1
    eeg = eeg1;
  else
    eeg = pop_mergeset(eeg, eeg1);
  end
end

% run ICA optional with our without PCA
% a window will pop-up as soon as ICA starts which allows to interrupt
% the ICA process. Please only press if you want to cancel the process
if cfg.PCA == 1
  eeg = pop_runica(eeg, 'icatype', 'runica', 'extended', 1, 'pca',...
    cfg.PCADIMS);
else
  eeg = pop_runica(eeg, 'icatype', 'runica', 'extended', 1);
end

% Run ICLabel
eeg = pop_iclabel(eeg,'Default');

% export ICA weights (to avoid duplicates of identical data)
ica = [];
ica.winv = eeg.icawinv;
ica.sphere = eeg.icasphere;
ica.weights = eeg.icaweights;
ica.chansind = eeg.icachansind;
ICLabel = eeg.etc.ic_classification.ICLabel;
save(fn_ica, 'ica','ICLabel')
diary(fn_log)

end
