function job = myspmeeg_import_ft(job)
% Imports FIELDTRIP data that is prep/epoched in MNE to SPM12
% 
% [USGAE]
% job = myspm_import_ft(job)
%
% [INPUT]
% job structure requires:
%  .fname_ft
%  .fname_rawmeg
%  .fname_out
%  .(overwrite)
%
% (cc) 2020, sgKIM.

if ~isfield(job,'overwrite')
  job.overwrite = false;
end

if exist(job.fname_out,'file') && ~job.overwrite
  ls(job.fname_out)
  return
end


% Import FIELDTRIP data to SPM:
load(job.fname_ft, 'epochs');
% epochs.grad.chanunit = ft_chanunit(epochs.grad, 'fT')
epochs.trial = single(epochs.trial)/1e-15; % converting T to fT
epochs.time = double(epochs.time);

[dn_out,~,~] = fileparts(job.fname_out);
[~,~] = mkdir(dn_out);
spm_eeg_ft2spm(epochs, job.fname_out);
trialinfo = epochs.trialinfo; % leave trial information for later use
clear epochs

% Load it to change headers:
load(job.fname_out, 'D') % saved in Matlab structure

% Correct MEG channel units to 'fT';
for i = 1:numel(D.channels)
  D.channels(i).units = 'fT';
end

% Correct trial labels:
for i = 1:numel(D.trials)
  D.trials(i).label = num2str(trialinfo(i));
end

% Import CTF header from the raw data to get fiducial points:
S = struct('dataset',job.fname_rawmeg, 'mode','header');
raw = spm_eeg_convert(S);
D.fiducials = raw.fiducials;

% Save it:
save(job.fname_out,'D')

end
