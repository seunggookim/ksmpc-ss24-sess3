%{
1 = Allegri, 2 = Olafur, 3 = Oblivion, 4 = Explosions, 5 = Vltava, 6 = Hoppipolla, 7 = Band of brothers
1 = Sad, 2 = Joy, 3 = Touched, 4 = Beauty, 5 = Connection, 6 = Warmth, 7 = Chills

%}
DnOsf = '/Volumes/APFS-2TB/extendedhome/MATLAB Drive/ksmpc-ss24-sess03/data/raw';
DnStim = '/Volumes/APFS-2TB/Keeper/MyLibrary/Work/Daniela/projects/musbhv/Vouskoski2022';


ScaleNames = ["sad","joy","touched","beauty","connection","warmth","chills"];
Tbl = readtable([DnOsf,'/ContinuousMusic_TimeSeries.csv']);
DnOut = './';
for iSong = 1:7
  for jScale = 1:7
    SubTbl = Tbl( (Tbl.Song == iSong) & (Tbl.Scale == jScale), :);
    Val = [];
    for t = 1:max(SubTbl.Sec)
      % rescale the value from [1,5] to [0,4] 
      % because "1" was the default ("no-response") value
      Val(t,1) = mean(SubTbl.Rating(SubTbl.Sec==t) - 1);
    end
    Sec = (1:max(SubTbl.Sec))';
    FnCsv = sprintf('%s/%s_song%i.csv', DnOut, ScaleNames(jScale), iSong);
    writetable(table(Sec, Val, VariableNames=['Sec',ScaleNames(jScale)]), FnCsv)
  end
end


Fnames = findfiles([DnStim,'/*flac']);
for iSong = 1:7
  [y, fs] = audioread(Fnames{iSong});
  env = db(abs(hilbert(mean(y,2))));
  [Env,Sec] = resample(env, (0:numel(env)-1)'/fs, 1);
  FnCsv = sprintf('%s/env_song%i.csv', DnOut, iSong);
  writetable(table(Sec, Env), FnCsv) 
end
