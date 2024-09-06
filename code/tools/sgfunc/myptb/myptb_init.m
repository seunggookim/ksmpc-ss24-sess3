function Job = mypbt_init(Job)
% Job = mypbt_init(Job)
% (cc) 2023, dr.seunggoo.kim@gmail.com

if not(exist('Job', 'var')), Job = []; end
Job = defaultjob(struct(Machine='mpieaeeg1', Debug=true, ...
  AudioMode=1, AudioSrateHz=44100, AudioReqLatencyClass=3, AudioNumChannels=2, ...
  TextFont='Monospace', TextSize=18, TextColor=[0 0 0], FillRect=round(255*0.8*[1 1 1]), Vspace=2 ), ...
  Job, mfilename);

%{
REFERENCE:

http://psychtoolbox.org/docs/PsychPortAudio-Open
‘mode’ Mode of operation. Defaults to 1 == sound playback only. Can be set to 2 == audio capture, or 3 for simultaneous
capture and playback of sound. Note however that mode 3 (full duplex) does not work reliably on all sound hardware. On
some hardware this mode may crash hard! There is also a special monitoring mode == 7, which only works for full duplex
devices when using the same number of input- and output-channels. This mode allows direct feedback of captured sounds
back to the speakers with minimal latency and without involvement of your script at all, however no sound can be
captured during this time and your code mostly doesn’t have any control over timing etc. You can also define a audio
device as a master device by adding the value 8 to mode. Master devices themselves are not directly used to playback or
capture sound. Instead one can create (multiple) slave devices that are attached to a master device. Each slave can be
controlled independently to playback or record sound through a subset of the channels of the master device. This
basically allows to virtualize a soundcard. See help for subfunction ‘OpenSlave’ for more info.

‘reqlatencyclass’ Allows to select how aggressive PsychPortAudio should be about minimizing sound latency and getting
good deterministic timing, i.e. how to trade off latency vs. system load and playing nicely with other sound
applications on the system. Level 0 means: Don’t care about latency or timing precision. This mode works always and with
all settings, plays nicely with other sound applications. Level 1 (the default) means: Try to get the lowest latency
that is possible under the constraint of reliable playback, freedom of choice for all parameters and interoperability
with other applications. Level 2 means: Take full control over the audio device, even if this causes other sound
applications to fail or shutdown. Level 3 means: As level 2, but request the most aggressive settings for the given
device. Level 4: Same as 3, but fail if device can’t meet the strictest requirements.

%}

switch Job.Machine
  case 'myubuntu22'
    setenv('WAYLAND_DISPLAY'); Screen('Preference','ConserveVRAM', 2^19);  % My lovely Ubuntu 22.04.1 uses XWayland
    Screen('Preference', 'SkipSyncTests', 1);  % Skip the test as XWayland is not officially supported
    AUDIO_DEVICE = 'sysdefault';

  case 'mpieaeeg1'
    %% MPI-EA LAB
    % REF: https://confluence.ae.mpg.de/display/LABWIK/Matlab
    %
    %[EEG1]
    % Presentation Computer  (PC-31):
    %  - Fujisu CELSIUS M740B, Xeon 3.5 GHz, 16 GB, SSD 1 TB
    %  - Windows 10 Pro (64-bit)
    %  - MATLAB R2021b (v9.10.0.1602886)
    %  - Psychphysics Toolbox 3 (v3.0.12)
    %  - external sound card: RME Fireface UCX
    %  - Loudspeakers: Neumann Model KH120A
    %
    %[Other equipments]
    % - Audiometer: MAICO MA25 (Neurosci department?)
    % - Sound level meter: RS PRO RS-1150 (ArtLab)

    %% TRIGGER OUTPUT
    %
    % LAB requires "io64.mexw64" and "inpoutx64.dll" to send triggers from PRESENTATION PC to EEG-RECORDING PC via
    % USB-ADAPTOR-BOX. These files are now in sgfunc/mypbt/external/
    %
    %```
    % Initialize LPT
    % ioObj = io64;
    % status = io64(ioObj);
    % address = hex2dec('dfff8');
    %
    % %Send trigger
    % triggervalue=255
    % io64(ioObj,address,triggervalue);
    % WaitSecs(0.005);
    % io64(ioObj,address,0);
    %```
    %
    % EEG1:LPT3:0xDFF8
    % EEG2:LPT3:0xFFF8
    % PSYPHY:LPT3:0xDFE8
    Job.PtrTrg = io64;
    Job.TrgStatus = io64(Job.PtrTrg);
    Job.TriggerAddress = hex2dec('dff8'); % EEG1


    %% INPUT DEVICES

    %% OUTPUT DEVICES

    % EXTERNAL SOUND CARD:
/mnt/projekte/2023-0360-MUSAFX/
    % DISPLAY:
    

  case 'mpieaeeg2'
    Job.TriggerAddress = hex2dec('fff8'); % EEG2
    Job.PtrTrg = io64;
    Job.TrgStatus = io64(Job.PtrTrg);

  case 'mpieapsyphy'
    Job.TriggerAddress = hex2dec('dfe8'); % PsyPhy
    Job.PtrTrg = io64;
    Job.TrgStatus = io64(Job.PtrTrg);


  case 'bic'

  case 'cobic'

  otherwise

end

%% VIDEO OUTPUT

% Set screen size: debug (640x320), non-debug (full screen, no cursor)
if Job.Debug
  ScreenRect = [0 0+50 640 320+50]; % ONLY FOR DEBUGGING
  Job.TextSize = 18;
else
  ScreenRect = []; % FULL SCREEN (should be mirrored to the projector)
  Job.TextSize = 24;
  %   HideCursor;
  %   cObj = onCleanup(@ShowCursor); % HAVE THIS GUY INSTEAD A FUNCTION (not the script...)
  %   IS RUNNING A SCRIPT SAFER BECAUSE YOU WILL STILL HAVE ALL LOCAL VARIABLES?
  %   WHY DON'T I CREATE A CLEANUP TASK TO SAVE EVERYTHING WHEN ABORTING? (A GOOD IDEA?)
end

% Create the pointer for the Window
Job.PtrWin = Screen('OpenWindow', 0, 0, ScreenRect);

% Set the Window properties
Screen('FillRect', Job.PtrWin, Job.FillRect);  % Screen background color
Screen('TextSize', Job.PtrWin, Job.TextSize);
Screen('TextFont', Job.PtrWin, Job.TextFont);


%% AUDIO OUTPUT

InitializePsychSound();  % Initialize the sound driver
PsychPortAudio('Close');  % close audio device
PsychPortAudio('Verbosity', 10);  % set the verbosity high
Dev = PsychPortAudio('GetDevices');  % get audio device IDs
Idx = find(contains(lower({Dev.DeviceName}), AUDIO_DEVICE));
if numel(Idx)==1
  logthis('Audio output device "%s" will be used...\n', Dev(Idx).DeviceName)
  Job.PtrAud = PsychPortAudio('Open', Dev(Idx).DeviceIndex, Job.AudioMode, Job.AudioReqLatencyClass, Job.AudioSrateHz, ...
    Job.AudioNumChannels);
else
  logthis('Which audio output device do know mean? AUDIO_DEVICE="%s", and the found devices are:\n', AUDIO_DEVICE)
  disp({Dev.DeviceName}')
  error('AUDIO DEVICE NOT DETEREMINED')
end


%% KEYBOARD/MOUSE INPUT
MaxPriority('GetSecs', 'KbWait', 'WaitSecs', 'KbEventAvail');






end