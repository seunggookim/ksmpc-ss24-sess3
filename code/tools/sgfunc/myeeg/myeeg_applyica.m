function eeg = myeeg_applyica(eeg, ica)
% eeg = myeeg_applyica(eeg, ica)
if ~isstruct(eeg)
  eeg = pop_loadset(eeg);
end
if ~isstruct(ica)
  load(ica,'ica');
end
eeg.icawinv = ica.winv;
eeg.icasphere = ica.sphere;
eeg.icaweights = ica.weights;
eeg.icachansind = ica.chansind;
eeg.setname = [eeg.setname '_ica'];
eeg = eeg_checkset(eeg);
end
