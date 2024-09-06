function EEG = myeeg_prep(CFG, EEG)

if ~isfield(CFG,'bpfreq'), CFG.bpfreq = [1 30]; end
% if ~isfield(CFG,'dsfreq'), CFG.dsfreq = 125; end
if ~exist('EEG','var')
  EEG = pop_loadset('filename',[CFG.fn_set,'.set'], 'filepath',CFG.dn_data);
end

%% Review possible SERIOUSLY (|Z|>10) bad channeles:
idx_badchan = [];
[~,idx,~] = pop_rejchan(EEG, 'elec',1:EEG.nbchan ,'threshold',10,...
  'norm','on', 'measure','kurt');
idx_badchan = [idx_badchan; idx];
[~,idx,~] = pop_rejchan(EEG, 'elec',1:EEG.nbchan ,'threshold',10,...
  'norm','on', 'measure','prob');
idx_badchan = [idx_badchan; idx];
[~,idx,~] = pop_rejchan(EEG, 'elec',1:EEG.nbchan ,'threshold',10,...
  'norm','on', 'measure','spec','freqrange',[.5 50]);
idx_badchan = unique([idx_badchan; idx]);
if numel(idx_badchan)
  warning('BAD CHANNELS FOUND: ');
  disp(cat(2,{EEG.chanlocs(idx_badchan).labels}))
  % now mark these channels, discard from PCA, and interpolate.
  return
end
CFG.prefix = CFG.fn_set;

%% Line noise suppression (but could be insufficient)
EEG = pop_cleanline(EEG, ...
  'LineFrequencies', [CFG.linefreq:CFG.linefreq:EEG.srate/2]);


%% BPF
EEG = pop_eegfiltnew(EEG, ...
  'lowcutoff', CFG.bpfreq(1), 'hicutoff', CFG.bpfreq(2));
EEG = eeg_checkset(EEG);
CFG.prefix = [CFG.prefix,'_',num2str(CFG.bpfreq(1)),'-',num2str(CFG.bpfreq(2)),'Hz'];
% pop_saveset(EEG,'filename',[CFG.prefix,'.set'], 'filepath',CFG.dn_data);


%% RUNICA on high-passed & downsamped data
TMP = EEG;
TMP = pop_eegfiltnew(TMP, 1, [], []); % high-pass > 1Hz
TMP = pop_resample(TMP, TMP.srate/8); % 1/8 downsampling
ICA = pop_runica(TMP, 'icatype','runica','dataset',1,...
  'options',{'extended' 1 'PCA' TMP.nbchan});
clear TMP


%% Apply ICA to the original data:
EEG.icaweights = ICA.icaweights;
EEG.icachansind = ICA.icachansind;
EEG.icasphere = ICA.icasphere;
EEG = eeg_checkset(EEG, 'ica');
CFG.prefix = [CFG.prefix,'_ica'];
pop_saveset(EEG,'filename',[CFG.fn_set,'_ica.set'], 'filepath',CFG.dn_data);


%% Component selection (ADJUST):
EEG = myeeg_autoadjust(fullfile(CFG.dn_data,[CFG.fn_set,'_ica.set']));


end