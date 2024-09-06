function Job = mypbt_test(Job)


%% DISPLAY TEST1
ScreenText = ['A puretone (1 kHz, 1 sec) is presented.\n',...
  '*PRESS ANY KEY IF YOU HEARD IT*'];
DrawFormattedText(Job.PtrWin, ScreenText, 'center', 'center', Job.TextColor, [], [], [], Job.Vspace);
Screen(Job.PtrWin,'Flip');

%% SOUND TEST: 1 sec 1kHz
Y = sin(linspace(0,1,Job.AudioSrateHz)*1000*pi)*0.2;
PsychPortAudio('Fillbuffer', Job.PtrAud, [Y; Y]);
PsychPortAudio('Volume', Job.PtrAud, 1);
PsychPortAudio('Start', Job.PtrAud, 1, 0, 0);

%% KEYBOARD TEST
KbReleaseWait;  % wait for all keys released
TARGET_KEYSTATUS = 0;  % 0=pressed, 1=released, 2=released>pressed, 3=released>pressed>released
Now_sec = GetSecs;
[Pressed_sec, IsKey] = KbWait([], TARGET_KEYSTATUS);

%% DISPLAY TEST2
ScreenText = sprintf('*KEY PRESSED*\nKeyId=[%i], Time=%f sec', find(IsKey), Pressed_sec-Now_sec);
DrawFormattedText(Job.PtrWin, ScreenText, 'center', 'center', Job.TextColor, [], [], [], Job.Vspace);
Screen(Job.PtrWin,'Flip');

%% MOUSE TRACE TEST
KbReleaseWait;  % wait for all keys released

DevIdx  = GetMouseIndices; % In Windows, all mouses are unified.
DevIdx = DevIdx(1);        % So, just pick the first one.

% "Create a keyboard queue for the first 3 valuators (x-axis, y-axis, scroll-wheel) in raw motion event data (4)"
% Okay but I need the desktop pixel coordinates...?????? 🤔
KbQueueCreate(DevIdx, [], 3, [], 4);

% PLAY MUSIC
KbQueueStart(DevIdx); % start queueing events

% init dT, x, y
SetMouse(0, 0);
[X, Y] = GetMouse;
oldTime = [];
dT = 0;

fid = fopen(FnameMouse, 'w');
fprintf(fid, 'dT[msecs]\txc\tyc\txi\tyi\tvx\tvy\twheel\n'); % WRITE the header line
fclose(fid); % flush

while ()
  while KbEventAvail(DevIdx)  % are there any events?
    Event = KbEventGet(d);  % READ the event
    if Event.Type == 1  % Is the event MOTION?

      % Accumulate absolute mouse position x,y from provided dx,dy movements:
      X = X + evt.Valuators(1);
      Y = Y + evt.Valuators(2);

      if ~isempty(oldTime)
        dT = Event.Time - oldTime;
      end
      oldTime = Event.Time;

      fid = fopen(FnameMouse, 'a'); % append a new line
      fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', dT*1000, Event.X, Event.Y, X, Y, Event.Valuators(1:3));
      fclose(fid);

    end
  end
end

% BAIL OUT: END of MUSIC or ESCAPE KEY PRESSED

KbQueueStop(DevIdx);  % stop queueing events
KbQueueRelease(DevIdx); % clean up the queue


end