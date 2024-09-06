function EEGout = myeeg_autoadjust (fname_in)
%
%
% Lines extracted from interface_ADJ
%
EEG = pop_loadset(fname_in);
[p1,f1,e1] = fileparts(fname_in);
fname_out = fullfile(p1,[f1,'_AUTOADJ',e1]);
fname_report = fullfile(p1,[f1,'_ADJUSTreport.txt']);
fname_list = fullfile(p1,[f1,'_ADJUSTlist.mat']);
givewarning
if exist(fname_out,'file')
  ls(fname_out)
  if nargout
    EEGout = pop_loadset(fname_out);
  end
  return
end

%% Run ADJUST

% Epoching

% ----------------------------------------------------
% | NOTE: epochs are extracted ONLY to make          |
% | ArtifactADJUST run                               |
% ----------------------------------------------------

% Check that ICA has been computed
if isempty(EEG.icaweights)
  warning('Please compute ICA before running ADJUST.');
  return;
end;

% Check that number of ICs = number of channels
if size(EEG.icaweights,1)<size(EEG.icaweights,2)
  disp(' ');
  warning('Number of ICs < number of channels.')
  disp(' ');
  disp('If ICA was not run on all channels, remove the excluded channels before running ADJUST.');
  disp('They can be reintroduced by interpolating them from neighbour channels after artifact removal.');
  disp('If some ICs were removed, please go back one step and run ADJUST on all ICs.')
  disp('The reduced IC statistics could lead to erroneous results.');
  return;
end;

if size(EEG.data,3)==1 % epochs must be extracted
  lag=5; %epoch duration (in seconds)
  disp(['Continuous dataset epoched in ' num2str(lag) ' sec long epochs to compute feature statistics over time']);
  % check whether '5sec' events are present
  si=0;
  for i=1:length(EEG.event)
    if(strcmp(EEG.event(1,i).type, '5sec')==1)
      si=1;
    end
  end
  if si==0 %add events
    ntrials=floor((EEG.xmax-EEG.xmin)/lag);
    nevents=length(EEG.event);
    for index=1:ntrials
      EEG.event(index+nevents).type='5sec';
      EEG.event(index+nevents).latency=1+(index-1)*lag*EEG.srate; %EEG.srate is the sampling frequency
      latency(index)=1+(index-1)*lag*EEG.srate;
    end;
    
    EEG=eeg_checkset(EEG,'eventconsistency');
  end
  
  EEGep = pop_epoch( EEG, {  '5sec'  }, [0 lag], 'newname', [EEG.setname '_ep5'] , 'epochinfo', 'yes');
  %         % removing baseline
  %         EEGep = pop_rmbase( EEGep, [0  0]);
  EEGep = eeg_checkset(EEGep);
  
  % collects ICA data from EEG
  if isempty(EEGep.icaact)
    disp('Warning: EEG.icaact missing! Recomputed from EEG.icaweights, EEG.icasphere and EEG.data');
    % Next instruction: see eeg_checkset
    EEGep.icaact = reshape(EEGep.icaweights*EEGep.icasphere*reshape(EEGep.data(1:size(EEGep.icaweights,1),:,:),[size(EEGep.icaweights,1) size(EEGep.data,2)*size(EEGep.data,3)]),[size(EEGep.icaweights,1) size(EEGep.data,2) size(EEGep.data,3)]);
  end;
  
  % Now that dataset is epoched, run ADJUST
  [art, horiz, vert, blink, disc,...
    soglia_DV, diff_var, soglia_K, med2_K, meanK, soglia_SED, med2_SED, SED, soglia_SAD, med2_SAD, SAD, ...
    soglia_GDSF, med2_GDSF, GDSF, soglia_V, med2_V, maxvar, soglia_D, maxdin]=ADJUST (EEGep, fname_report);
else % data are epoched... let's work on the existing epochs
  
  % collects ICA data from EEG
  if isempty(EEG.icaact)
    disp('Warning: EEG.icaact missing! Recomputed from EEG.icaweights, EEG.icasphere and EEG.data');
    % Next instruction: see eeg_checkset
    EEG.icaact = reshape(EEG.icaweights*EEG.icasphere*reshape(EEG.data(1:size(EEG.icaweights,1),:,:),[size(EEG.icaweights,1) size(EEG.data,2)*size(EEG.data,3)]),[size(EEG.icaweights,1) size(EEG.data,2) size(EEG.data,3)]);
  end;
  
  % run ADJUST
  [art, horiz, vert, blink, disc,...
    soglia_DV, diff_var, soglia_K, med2_K, meanK, soglia_SED, med2_SED, SED, soglia_SAD, med2_SAD, SAD, ...
    soglia_GDSF, med2_GDSF, GDSF, soglia_V, med2_V, maxvar, soglia_D, maxdin]=ADJUST (EEG, fname_report);
end
%% Saving artifacted ICs for further analysis
eb=blink;
hem=horiz;
vem=vert;
gd=disc;
aic=unique([blink disc horiz vert]);
save (fname_list, 'eb', 'hem', 'vem', 'gd' , 'aic');
disp(['Artifact ICs list saved in ' fname_list]);
disp(' ')

%% Creating component plots
fname_png = fullfile(p1,[f1,'_ADJ.png']);
if ~exist(fname_png, 'file')
  figure('position',[1921, 78, 1470, 896], 'visible','off')
  myeeg_plotcomp(EEG, aic)
  export_fig(fname_png)
end

%%
disp(['Removing component(s): ',num2str(aic)])
EEGout = pop_subcomp(EEG, aic);
pop_saveset(EEGout, fname_out);
fprintf('Saving in ')
ls(fname_out)
if ~nargout
  clear EEGout
end
givewarning
end

function givewarning
warning('***THIS IS ONLY FOR PRELIMINARY ANALYSIS!!!!***')
warning('***YOU HAVE TO INSPECT ALL COMPONENTS CAREFULLY FOR ACTUAL ANALYSIS!!!***')
end
