function eeg = myeeg_loadwithica(fname)
% eeg = myeeg_loadwithica(fname)
% (cc) 2012, sgKIM

[p1,f1,~] = fileparts(fname);
fn_ica = fullfile(p1,[f1,'_ica.mat']);
if ~isfile(fn_ica)
  error('ICA FILE (%s) NOT FOUND.', fn_ica)
end
eeg = pop_loadset(fname);
load(fn_ica,'ica');
eeg.icawinv = ica.winv;
eeg.icasphere = ica.sphere;
eeg.icaweights = ica.weights;
eeg.icachansind = ica.chansind;
eeg.setname = [eeg.setname '_ica'];
eeg = eeg_checkset(eeg);
end
