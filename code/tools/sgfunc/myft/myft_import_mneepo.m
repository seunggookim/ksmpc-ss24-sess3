function epochs = myft_import_mneepo (fn_epo, fn_rawmeg, dn_out, overwrite)
% converts MNE epoched structure to FieldTrip timelocked structure.
%
% [USAGE]
% epochs = myft_import_mneepo (fn_epo, fn_rawmeg, dn_out, [overwrite])
%
% [INPUT]
% fn_epo      '1xN' filename for an epoched MEG file processed by MNE
% fn_rawmeg   '1xN' filename for an original MEG raw data
% dn_out      '1xN' dirname to save fieldtrip format .MAT file
% (overwrite) [1x1] logical flag to overwrite files (default=false)
%
% [OUTPUT]
% epochs      (1x1) imported FieldTrip structure of epoched MEG data
%
% (cc) 2020, sgKIM, solleo@gmail.com

%{
NOTE: Basically ft_preprocessing can read MNE epoched fiff files but it is
extremely slow (takes 10 min to read 720 epochs). So we'll just use it to
read header and replace the TRIAL field.
%}

if ~exist('overwrite','var')
  overwrite = false;
end
[~,f1,~] = myfileparts(fn_epo);
[~,~] = mkdir(dn_out);
fn_out = fullfile(dn_out,[f1,'.mat']);
if exist(fn_out,'file') && ~overwrite
  ls(fn_out)
  return
end

%% SENSOR INFO
% FN_RAWMEG should be run1, to which all runs are coregistered.
grad = ft_read_sens(fn_rawmeg, 'senstype','meg');
grad = ft_convert_units(grad, 'cm');
elec = ft_read_sens(fn_rawmeg, 'senstype','eeg');
elec = ft_convert_units(elec, 'cm');


%% EVENT INFO

% Get coordinate system info from the raw file:
%raw = ft_preprocessing(struct('dataset',fn_raw, 'trl',[1 1 0]));

% Read event info & data from the epoched file:
epo = ft_read_header(fn_epo);%, 'coordsys', raw.grad.coordsys);

% Parse event IDs:
str = strsplit(epo.orig.epochs.event_id,';');
eventid = cell(1,numel(str));
for i = 1:numel(str)
  temp = strsplit(str{i},':');
  eventid{str2double(temp{2})} = temp{1};
end

% Trial definition:
trl = double(epo.orig.epochs.events);
trl(:,2) = trl(:,1) + epo.nSamples - 1;


%% DATA and HEADER (tSSS'ed) INFO

% Import epoched data (only one sample; somehow it's extremely slow)
epochs = ft_preprocessing(struct('dataset',fn_epo,'trl',[1 1 0]));
epochs.time = cell(1,epochs.hdr.nTrials);
epochs.trial = cell(1,epochs.hdr.nTrials);
for i = 1:epochs.hdr.nTrials
  epochs.time{i} = single(epo.orig.epochs.times);
  epochs.trial{i} = single(squeeze(epo.orig.epochs.data(i,:,:)));
end
epochs.sampleinfo = trl(:,1:2);
epochs.trialinfo = trl(:,3);
epochs.cfg.eventid = eventid;
epochs.cfg.mneinfo = epo.orig.epochs.info;
epochs.grad = grad;
epochs.elec = elec;

warning('Overwriting channel labels with grad.label')
epochs.label = grad.label(1:numel(epochs.label));

% Validate FT format:
epochs = ft_checkdata(epochs, 'feedback','yes', 'datatype','timelock');
epochs.trial = single(epochs.trial);

% Save data:
save(fn_out,'epochs')
fprintf('Saving ')
ls(fn_out)
if ~nargout, clear epochs; end


end