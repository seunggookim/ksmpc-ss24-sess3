function eeg = myeeg_iclabel(eeg, thres)
% eeg = myeeg_iclabel(fname)
% eeg = myeeg_iclabel(eeg)
% eeg = myeeg_iclabel(fname, thres)
% eeg = myeeg_iclabel(eeg, thres)

if ~isstruct(eeg)
  eeg = myeeg_loadwithica(eeg);
end
if ~exist('thres','var')
  thres = 0.8;
end
eeg = pop_iclabel(eeg, 'default');
ICLabel = eeg.etc.ic_classification.ICLabel;
[~,f1,e1] = fileparts(eeg.filename);
save(fullfile(eeg.filepath, [f1,'_ICLabel.mat']), 'ICLabel')

eeg = pop_icflag(eeg, [NaN NaN;thres 1;thres 1;thres 1;thres 1;thres 1;NaN NaN]);
eeg = pop_subcomp(eeg, []);
eeg.setname = [eeg.setname,'_proj'];
eeg.filename = [f1,'_proj',e1];
pop_saveset(eeg, 'filename',eeg.filename, 'filepath',eeg.filepath, ...
  'savemode','onefile');

% SEE `eeg.etc.ic_classification` for classification results
end