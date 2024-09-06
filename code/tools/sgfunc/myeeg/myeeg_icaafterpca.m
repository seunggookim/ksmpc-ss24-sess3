function ica = myeeg_icaafterpca(fn_run, cfg)
% ica = myeeg_icaafterpca(fn_run, cfg)

if ~exist('cfg','var')
  cfg =[];
  cfg.HP = 1;         % cut-off frequency high-pass filter [Hz] only for ICA
  cfg.LP = 40;        % cut-off frequency low-pass filter [Hz] only for ICA
  cfg.SRATE = 250;    % downsample data for ICA
  cfg.HP_ord = 500;   % high-pass filter order depends on sampling rate
  cfg.LP_ord = 100;   % low-pass filter order depends on sampling rate
  cfg.PRUNE = 3;      % artifact rejection threshold in SD (for ICA only)
  cfg.PCA = 1;        % choose PCA option for ICA
  cfg.PCADIMS = 50;   % PCA dimension if PCA option is true
end

% load rawdata (already saved as set file)
eeg = pop_loadset(fn_run);
[p1,f1,e1] = fileparts(fn_run);

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
eeg = pop_jointprob(eeg, 1, [1:size(eeg.data,1)], ...
  cfg.PRUNE, cfg.PRUNE, 0, 1, 0);

% run ICA optional with our without PCA
% a window will pop-up as soon as ICA starts which allows to interrupt
% the ICA process. Please only press if you want to cancel the process
if cfg.PCA == 1
  eeg = pop_runica(eeg, 'icatype', 'runica', 'extended', 1, 'pca',...
    cfg.PCADIMS);
else
  eeg = pop_runica(eeg, 'icatype', 'runica', 'extended', 1);
end

% export ICA weights (to avoid duplicates of identical data)
ica = [];
ica.winv = eeg.icawinv;
ica.sphere = eeg.icasphere;
ica.weights = eeg.icaweights;
ica.chansind = eeg.icachansind;
save(fullfile(p1,[f1,'_ica.mat']),'ica')

end
