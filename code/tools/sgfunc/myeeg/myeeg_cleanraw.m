function EEG = myeeg_cleanraw(FnameSet)
% Eeg = myeeg_cleanraw(FnameSet)
if isstruct(FnameSet)
  EEG = FnameSet;
else
  EEG = pop_loadset(FnameSet);
end
chanlocs = EEG.chanlocs;

% raw-clean (correct bursts):
EEG = pop_clean_rawdata(EEG, 'Highpass',[.25 .75], 'WindowCriterion', 0.25, ...
  'FlatlineCriterion',5, 'LineNoiseCriterion',4, 'ChannelCriterion',0.8, 'BurstCriterion',20, 'BurstRejection','off');

% interpolate removed channels:
if not(isempty([chanlocs.X])) && not(isempty(setdiff({chanlocs.labels}, {EEG.chanlocs.labels})))
  EEG = pop_interp(EEG, chanlocs, 'spherical');
end

% % save:
% EEG.setname = [EEG.setname,'_cleaned'];
% [~,f1,e1] = fileparts(EEG.filename);
% EEG.filename = [f1,'_cleaned',e1];
% pop_saveset(EEG, 'filename',EEG.filename, 'filepath',EEG.filepath, 'savemode','onefile');

end