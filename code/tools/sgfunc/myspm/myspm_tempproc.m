function Job = myspm_tempproc(Job)


%% Read parameters from JSON, DICOM, NIFTI headers or user inputs:
if not(isempty(Job.FnameJson))
  Job.FnameJson = [Job.EpiDname,'/',Job.EpiPrefix,'.json'];
  ls(Job.FnameJson)
end

% [path1,Name1,ext1] = myfileparts(Job.FnameEpi);
% FnameEpi = [path1,'/',Name1,ext1];
% ls(FnameEpi)

%% find TR
if isempty(Job.Tr_sec) % manually given?
  if isfile(Job.FnameJson) % using json files
    logthis('Found json file: '); ls(Job.FnameJson);
    logthis('Reading TR & SliceTiming from json\n');
    Json = jsondecode(fileread(Job.FnameJson));
    Job.Tr_sec = Json.RepetitionTime; % *** it's SEC!! ***
  elseif not(isempty(Job.FnameDcm))
    if isfile(fn_dcm)
      logthis('Found dicom file: '); ls(Job.FnameDcm);
      Hdr2 = spm_dicom_headers(Job.FnameDcm);
      Job.Tr_sec  = Hdr2{1}.RepetitionTime/1000; % it's msec!!!
    end
  else
    error('TR cannot be determined!')
  end
end
logthis('TR = %i msec\n',Job.Tr_sec*1000);

%% CROP dummy scans (for old scanners or new protocols)
if not(isempty(Job.Dummy_sec))
  FnameEpi_old = Job.FnameEpi;
  n = ceil(Job.Dummy_sec / Job.Tr_sec);
  fprintf('Removing first %i sec (%i volumes)\n',Job.Dummy_sec, n);
  FnameEpi_new = [Job.EpiDname,'/',Job.EpiPrefix,'_skip',num2str(n),ext1];
  if ~exist(FnameEpi_new,'file')
    V = spm_vol(FnameEpi_old);
    for i = (n+1):numel(V)
      Vi = V(i);
      Y = spm_read_vols(Vi);
      Vi.fname = FnameEpi_new;
      Vi.n(1) = i-n;
      spm_write_vol(Vi,Y);
    end
  end
  ls(FnameEpi_new)
  Job.FnameEpi = FnameEpi_new;
  [~,Job.EpiPrefix,~] = fileparts(Job.FnameEpi);
end

FnameOut = fullfile(Job.EpiDname, ['a',Job.EpiPrefix,'.nii']);
if isfile(FnameOut)
  logthis('File found: '); ls(FnameOut)
  return
end

%% FIND slice timing
V = spm_vol(Job.FnameEpi);
NumFrames = numel(V);
logthis('Number of frames = %i\n',NumFrames);

if isfile(Job.FnameJson)
  if contains(upper(Json.Manufacturer),'SIEMENS')
    logthis('Siemens scanner detected.\n')
    slice_order_msec = Json.SliceTiming*1000;
    ref_slice_msec = Job.Tr_sec*1000/2;
    if ~isfield(Json,'TotalReadoutTime') && isfield(Json,'EffectiveEchoSpacing')
      Json.TotalReadoutTime = 1/(1/Json.EffectiveEchoSpacing/Json.BaseResolution);
      % REF: https://lcni.uoregon.edu/kb-articles/kb-0003 (SPM definition)
    end
    totalreadout_msec =  Json.TotalReadoutTime * 1000;
  elseif contains(Json.Manufacturer,'GE')
    logthis('GE scanner detected.\n')
    if isfield(Job,'FnameBxh') % for duke-BIAC
      logthis('GE-BIAC?\n')
      bxh = parsebxh(Job.FnameBxh);
      slice_order_msec = bxh.acquisitiontime_msec;
      ref_slice_msec = Job.Tr_sec*1000/2;
    else % other GE scanners:
      slice_order_msec = Json.SliceTiming;
      ref_slice_msec = Job.Tr_sec*1000/2;
      totalreadout_msec = Json.TotalReadoutTime * 1000;
    end
  end
else
  % if no json
  if Job.Tr_sec >= 6 && ~isfield(Job,'noSTC') && ~isfield(Job,'FnameDcm')
    warning(...
      ['TR is long (>= 6 s) and found no dicom file (.FnameDcm)',...
      'to read actual slice time. Thus setting .noSTC=1']);
    logthis(['[!] It may possible to use vectors in .slice_order_msec ',...
      'and .ref_slice_msec\n']);
    Job.noSTC=1;
  end
  % find slice timing in msec and repetition time in sec from a example
  % DICOM
  if isfield(Job,'FnameDcm')
    Hdr2 = spm_dicom_headers(Job.FnameDcm);
    if isfield(Hdr2{1},'Private_0019_1029') % recent Siemens scanners
      slice_order_msec = Hdr2{1}.Private_0019_1029; % in msec
    else
      if ~isfield(Job,'slice_order')
        error('No slice timing information in the DICOM header!')
      end
    end
    ref_slice_msec = Job.Tr_sec*1000/2;
  elseif isfield(Job,'slice_order_msec')
    disp('Applying manually entered .slice_order_msec');
    slice_order_msec = Job.slice_order_msec;
    ref_slice_msec = Job.Tr_sec*1000/2;
    Job.noSTC = 0;
  else
    warning(['NO slice timing information is given. Skipping STC ', ...
      'CREATING a link...']);
    Job.noSTC = 1;
  end
end
% override of manually given:
if isfield(Job,'slice_order_msec')
  slice_order_msec = Job.slice_order_msec;
end
if isfield(Job,'ref_slice_msec')
  ref_slice_msec = Job.ref_slice_msec;
end


if isfield(Job,'noSTC') && Job.noSTC
  logthis(['Create a link for ',Job.FnameEpi,' as ',FnameOut,'..\n']);
  unix(['ln -sf ',Job.FnameEpi,' ',FnameOut])
else
  logthis('Slice order (msec) =\n'); 
  disp(reshape(slice_order_msec,1,[]));
  logthis('Reference slice (msec) = %.4f\n',ref_slice_msec);
  
  if Job.UseFslStc
    if isfield(Job,'FnameDcm')
      Hdr2 = spm_dicom_headers(Job.FnameDcm);
      slice_order_msec = Hdr2{1}.Private_0019_1029; % in msec
      T = slice_order_msec'./(Job.Tr_sec*1000); % in TR
    else
      T = (Job.slice_order_msec-1)./max(Job.slice_order_msec);
    end
    fname_stc = [p_epi,'/',f1,'_timing.txt'];
    dlmwrite(fname_stc,T);
    unix(['FSLOUTPUTTYPE=NIFTI; slicetimer.fsl -i ',Job.FnameEpi,...
      ' -o ',FnameOut,...
      ' --tcustom=',fname_stc,' -r ',num2str(Job.Tr_sec),' -v -d 3'])
  else
    Job1 = struct('fname_epi',Job.FnameEpi, 'NumFrames',NumFrames, ...
      'TR_sec',Job.Tr_sec, 'slice_order_msec',slice_order_msec, ...
      'ref_slice_msec',ref_slice_msec);
    myspm_stc(Job1);
  end
  
end

end