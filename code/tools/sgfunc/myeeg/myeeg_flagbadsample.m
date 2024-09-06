function eeg = myeeg_flagbadsample(fn)
% EEG = myeeg_flagbadsample(fn)
if isstruct(fn)
  eeg = fn;
else
  eeg = pop_loadset(fn);
end
% detect BURST after low-pass filtering?
cleaned = pop_eegfiltnew(eeg, 1, 30);
cleaned = clean_artifacts(cleaned, ...
  'Highpass',[.25 .75], 'WindowCriterion', 0.25, ...
  'FlatlineCriterion',5, 'LineNoiseCriterion',4,  ...
  'BurstCriterion',20, 'BurstRejection','on');

eeg.etc.clean_sample_mask = cleaned.etc.clean_sample_mask;
fprintf('[%s] %i samples (%.2f%%) are marked as BAD\n',...
  [eeg.filepath filesep eeg.filename], sum(eeg.etc.clean_sample_mask),...
  100*mean(eeg.etc.clean_sample_mask));
pop_saveset(eeg,'savemode','resave')

end