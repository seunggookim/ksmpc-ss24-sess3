function [EEG, Cfg] = myeeg_importbv(Cfg)
% EEG = myeeg_importbv(Cfg)
%
% Cfg requires:
%  .FnameRaw        '1xN'  filename of the raw data (with any extension)
%  .FnameSet        '1xN'  filename to save in EEGLAB format (with any extension)
% (.SamplesToRead)  [1x2]  samples to read [onset, offset]
% (.ChannelsToRead) [1xN]  channels to read [chanIndex1, chanIndex2, ... chanIndexN]
% (.IsSaveset)      [1x1]  true (default) | false
%
% For Brain Vision EEG system. The data files are:
% ${FnameRaw}.eeg     Actual binary data
% ${FnameRaw}.vhdr    ASCII text header file (data format, channel info, ect)
% ${FnameRaw}.vmrk    ASCII text file for markers of stimuli and responses
%
% To import data into eeglab, first you need to install a plug-in that allows to you import Brain Vision data. From the
% EEGLAB menu:
%
% File > Manage EEGLAB Extensions > Data Import Extensions
%
% Then from a popup window, click on Install tick for "bva-io" which "Imports Brain Vision Analyser data files" then
% click OK to install.
%
% (cc) 2019-2024, dr.seunggoo.kim@gmail.com

if not(exist('eeg_checkset','file'))
  error('Initialize eeglab before running %s!', mfilename)
end
if not(exist('pop_loadbv','file'))
  error('Brain Vision data importing plug-in needs to be installed!')
  help(mfilename)
end
if isfile(Cfg.FnameSet)
  logthis('Already here: ')
  ls(Cfg.FnameSet)
  if nargout
    EEG = pop_loadset(Cfg.FnameSet);
  end
  return
end

logthis('File to import: ')
ls(Cfg.FnameRaw)
[Cfg.DnameRaw, Cfg.PrefixRaw, ~] = fileparts(Cfg.FnameRaw);
[Cfg.DnameSet, Cfg.PrefixSet, ~] = fileparts(Cfg.FnameSet);
Cfg = defaultcfg( struct(SamplesToRead=[], ChannelsToRead=[], IsSaveset=true, FnameChanlocs=[]), Cfg, mfilename);

%% import Brain Vision files
EEG = pop_loadbv(Cfg.DnameRaw, [Cfg.PrefixRaw,'.vhdr'], Cfg.SamplesToRead, Cfg.ChannelsToRead);
if not(isempty(Cfg.FnameChanlocs))
  EEG = pop_editset(EEG, 'chanlocs', Cfg.FnameChanlocs); % set channel locations
end
EEG = eeg_checkset(EEG); % check if data structure is valid

%% save in EEGLAB format
if Cfg.IsSaveset
  [~,~] = mkdir(Cfg.DnameSet);
  EEG = pop_saveset(EEG, 'filename',[Cfg.PrefixSet,'.set'], 'filepath',Cfg.DnameSet);
  logthis('Imported file: ')
  ls(Cfg.FnameSet)
end

end
