function [EEG, cfg] = myeeg_importdata_olab(cfg)
% EEG = myeeg_importdata(cfg) - imports O-lab's EEG data
%
% cfg requires:
%  .fn_raw     '1xN'  filename of the raw data (with any extension)
%  .fn_out     '1xN'  filename to save in EEGLAB format (with any extension)
% (.smpl2read) [1x2]  samples to read
% (.chan2read) [1x2]  channels to read
% (.saveset)   [1x1]  true (default) | false
%
% meta info...?
% (.subject)   '1xN'
% (.group)     '1xN'
% (.condition) '1xN'
% (.session)   [1x1]
%
% O-lab uses Brain Vision EEG system. The data files are:
% ${fn_raw}.eeg     Actual binary data
% ${fn_raw}.vhdr    ASCII text header file (data format, channel info, ect)
% ${fn_raw}.vmrk    ASCII text file for markers of stimuli and responses
%
% To import data into eeglab, first you need to install a plug-in that allows
% to you import Brain Vision data. From the EEGLAB menu:
%
% File > Manage EEGLAB Extensions > Data Import Extensions
%
% Then from a popup window, click on Install tick for
% "bva-io" which "Imports Brain Vision Analyser data files"
% then click OK to install.
%
% (cc) 2019, sgKIM. solleo@gmail.com

if ~isfield(cfg,'smpl2read'), cfg.smpl2read = []; end
if ~isfield(cfg,'chan2read'), cfg.chan2read = []; end
if ~exist('eeg_checkset','file')
  error('Initialize eeglab before running %s!', mfilename)
end
if ~exist('pop_loadbv','file')
  error('Brain Vision data importing plug-in needs to be installed!')
  help(mfilename)
end
if ~isfield(cfg,'dn_data')
  [cfg.dn_data,~,~] = fileparts(cfg.fn_raw);
end
[~,cfg.fn_raw,~] = fileparts(cfg.fn_raw);
%% import Brain Vision files
EEG = pop_loadbv(cfg.dn_data, [cfg.fn_raw,'.vhdr'], cfg.smpl2read, cfg.chan2read);  
EEG = eeg_checkset(EEG); % check if data structure is valid

%% set channel locations
[mypath,~,~] = fileparts(mfilename('fullpath'));
chanlabels = {EEG.chanlocs.labels};
num_eegchan = sum(contains(chanlabels,'Ch'));
num_auxchan = sum(contains(chanlabels,'Aux'));
if num_eegchan == 63 && num_auxchan == 0
  fn_chanloc = fullfile(mypath,'mw_RRNorth_channel_locs_63chans_ref32.ced');
elseif num_eegchan == 63 && num_auxchan == 1
  fn_chanloc = fullfile(mypath,'mw_RRNorth_channel_locs_63chans_ref32_1aux.ced');
elseif num_eegchan == 64 && num_auxchan == 0
  fn_chanloc = fullfile(mypath,'mw_RRNorth_channel_locs_64chans.ced');
else
  disp(chanlabels)
  warning('Channel locations cannot be automatically determined.');
end
if exist('fn_chanloc','var')
  EEG = pop_editset(EEG, 'chanlocs', fn_chanloc);
end
%{
Channel types: eeg, eog, ref, aux

Channel labels:
  LHEOG - left horizontal EOG
  LM    - left mastoid
  RM    - right mastoid (reference)

Channel location coordinate system: ALS

%}

%% meta info..?
if isfield(cfg,'group')
  EEG.group = cfg.group;
end
if isfield(cfg,'subject')
  EEG.subject = cfg.subject;
end
if isfield(cfg,'session')
  EEG.session = cfg.session;
end
if isfield(cfg,'condtion')
  EEG.condtion = cfg.condtion;
end

%% save in EEGLAB format
if ~isfield(cfg,'dn_out') && ~isfield(cfg,'fn_out')
  cfg.dn_out = cfg.dn_data;
elseif ~isfield(cfg,'dn_out') && isfield(cfg,'fn_out')
  [cfg.dn_out,~,~] = fileparts(cfg.fn_out);
end
[~,~] = mkdir(cfg.dn_out);
if ~isfield(cfg,'fn_out')
  [~,cfg.fn_out_prefix,~] = fileparts(cfg.fn_raw);
else
  [~,cfg.prefix_out,~] = fileparts(cfg.fn_out);
end
if ~isfield(cfg,'saveset')
  cfg.saveset = true;
end
if cfg.saveset
  EEG = pop_saveset(EEG, 'filename', [cfg.prefix_out,'.set'], ...
    'filepath',cfg.dn_out);
end
if ~nargout
  clear EEG cfg
end
end
