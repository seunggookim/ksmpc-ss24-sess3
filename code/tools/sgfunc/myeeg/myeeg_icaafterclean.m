function ica = myeeg_icaafterclean(fn_run, cfg)
% ica = myeeg_icaafterpca(fn_run, cfg)


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

% load rawdata (already saved as set file)
if isstruct(fn_run)
  eeg = fn_run;
  fn_run = fullfile(eeg.filepath, eeg.filename);
else
  eeg = pop_loadset(fn_run);
end
[p1,f1,~] = fileparts(fn_run);
fn_ica = fullfile(p1,[f1,'_ica.mat']);
if isfile(fn_ica)
  load(fn_ica, 'ica')
  return
end
chanlocs = eeg.chanlocs;

% Aggressively clean the raw data (only for ICA)
eeg = pop_clean_rawdata(eeg, ...
  'FlatlineCriterion',5, 'ChannelCriterion',0.8, 'LineNoiseCriterion',4,...
  'Highpass',[.25 .75], 'BurstCriterion',20, 'WindowCriterion',.25, ...
  'BurstRejection','on', 'Distance','Euclidian',...
  'WindowCriterionTolerances',[-Inf 7]);
eeg = pop_interp(eeg, chanlocs, 'spherical');

% apply low pass filter (only for ICA)
eeg = pop_eegfiltnew(eeg, 'hicutoff',cfg.LP, 'filtorder', cfg.LP_ord);
% downsample data - may be an optional step
eeg = pop_resample(eeg, cfg.SRATE);
% apply high pass filter (only for ICA)
eeg = pop_eegfiltnew(eeg, 'locutoff',cfg.HP, 'filtorder', cfg.HP_ord);

% create dummy events and epoch data to these dummy events (!)
eeg = eeg_regepochs(eeg, 'recurrence', 1, 'eventtype', '999');
eeg = eeg_checkset(eeg, 'eventconsistency');

% remove epochs with artefacts to improve ICA training
eeg = pop_jointprob(eeg, 1, [1:size(eeg.data,1)], cfg.PRUNE, cfg.PRUNE, 0, 1, 0);

% run ICA optional with our without PCA
% a window will pop-up as soon as ICA starts which allows to interrupt
% the ICA process. Please only press if you want to cancel the process
if cfg.PCA == 1
  eeg = pop_runica(eeg, 'icatype', 'runica', 'extended', 1, 'pca',cfg.PCADIMS);
else
  eeg = pop_runica(eeg, 'icatype', 'runica', 'extended', 1);
end

% export ICA weights (to avoid duplicates of identical data)
ica = [];
ica.winv = eeg.icawinv;
ica.sphere = eeg.icasphere;
ica.weights = eeg.icaweights;
ica.chansind = eeg.icachansind;
save(fn_ica,'ica')

end
