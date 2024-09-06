function Job = myspm_fmriprep12 (Job)
% job = myspm_fmriprep12 (job)
%
% For preprocessing of EPI and T1w images. This script is based on a batch
% script included in SPM12.
% It does slice-timing-correction before realignment/unwarping.
%
% CURRENTLY, THIS SCRIPT REQUIRES FSL AND FREESURFER TOO FOR TINY BUT
% VERY CONVENIENT UTILITIES.
%
% This does:
%  (1) STRUCTURAL DATA (FnameT1w)
%      -> unified segmentation+normalization -> brain-masking
%      -> normalize into mni
%  (2) FIELDMAP DATA (fname_pha/fname_mag)
%      -> construct VDM
%  (3) FUNCTIONAL DATA (FnameEpi)
%      -> slice timing correction -> unified unwarp+realign
%      -> coregister into t1w (only header)
%      -> normalize into mni (resampling) -> smoothing (if defined)
%
% Job structure (1x1) requires:
%:Data:
%  .FnameEpi   : a string for a single session
% (.FnameJson) : (or with the same filename as the .nii file to be found)
%  .FnameT1w   : a string for anatomical image (full path)
%
%:T1w segmentation:
% (.Segment)   : 'spm12' (default) | 'cat12' 
%
%:Registration:
% (.UseAnts)    : 0 [default] | 1 (takes LOOONG: need MYANTS wrappers)
%
%:Cropping:
% (.Dummy_sec)  : removing first K sec for unsteady state (adding suffix in
%                 TR)
%
%:Slice timing correction:
% * Tr_sec, slice_oder_msec, ref_slices_msec can be read from a paired json
% file of FnameEpi
%
% (.SliceOrder_msec) : can be 1-based indices or actual time from the
%                       pulse (msec)
% (.RefSlice_msec)   : reference slice (in the same unit as
%                       .slice_order_msec)
% (.FanmeDcm)  : DICOM file of a volume to read slice time (only for
%                 SIEMENS scanners)
% (.UseFslStc) : use FSL instead of SPM for a very big file that SPM can't
%                 handle! (But are you sure you're doing right with the
%                 nortorious FSL slice-timing syntax?)
% (.NoStc)      : for a sparse sampling (e.g. TR>TA) STC would make no
%                 effect or could introduce a strange interpolation between
%                 distant time points?
% (.Tr_sec)     : can be read from a paired .json file.
%
%:Fieldmap: geometric distortion correction:
% (.FnameMag)  : magnitude a short TE
% (.FnamePha)  : phase difference
% * TEs_fmap and totalreadout_msec can be read from paired json files of
% fname_mag and fname_pha
% (.TeFmap)
% (.TotalReadout_msec)
% (.Fnamevdm)  : if it's already computed
%
% [Unwarping+realignment]
% (always run with default parameters)
%
%:CompCor regressors:
%  .Bpf
%  .nPcs
%
%:Spatial normalization:
% (.Vox_mm)     : resampling resolution for EPI in MNI152 space
% (.Bbox_mm)
%
%:Smoothing:
% (.Fwhm_mm)    : smoothing in MNI152 space (default = none)
%
% This script also uses original and modified MATLAB code by others:
% - DPARSF by Chao-Gan YAN: http://rfmri.org/DPARSF
%
% (cc) 2015-2023. dr.seunggoo.kim@gmail.com

if nargin==0, help(mfilename); return; end

%% ---CHECK INPUTS---
assert(isfield(Job,'FnameEpi'))
Job.FnameEpi = getfullpath(Job.FnameEpi);
assert(isfile(Job.FnameEpi))

assert(isfield(Job,'FnameT1w'))
Job.FnameT1w = getfullpath(Job.FnameT1w);
assert(isfile(Job.FnameT1w))

[Job.EpiDname, Job.EpiPrefix, Ext] = fileparts(Job.FnameEpi);
assert(strcmpi(Ext,'.nii'), 'FnameEpi="%s" is not .nii!', Job.FnameEpi)

[Job.T1wDname, Job.T1wPrefix, Ext] = fileparts(Job.FnameT1w);
assert(strcmpi(Ext,'.nii'), 'FnameT1w="%s" is not .nii!', Job.FnameT1w)


%% ---Set default options---
V = spm_vol(Job.FnameEpi);
DefaultJob = struct('UseAnts',false, 'Fwhm_mm',[], 'Bpf',[1/128 inf], ...
  'nPcs',16, 'Bbox_mm',[], 'Dummy_sec',[], 'UseFslStc',false, ...
  'NoStc',false, 'FnameJson',[], 'FanmeDcm',[], 'Tr_sec',[], ...
  'Segment','spm12', ...
  'Vox_mm',[1 1 1]*round(abs(diag(V(1).mat(1:3,1:3)))));
Job = defaultjob(DefaultJob, Job, mfilename);


%% KEEP DIARY
Job.FnameLog = fullfile(Job.EpiDname, [mfilename,'_',Job.EpiPrefix,'.log']);
tStart = tic; logfile(Job.FnameLog, Job, true);


%% ---CHECK ENVIRONMENT---
%% 0-1. Windows
if ispc
  error('Use Linux or Mac (I don''t have a Windows machine to test)')
end

%% 0-2. FSL
assert(isfile(fullfile(getenv('FSLDIR'),'bin','fslmaths')), ...
  'FSL not accessible!')

%% 0-3. Freesurfer
assert(isfile(fullfile(getenv('FREESURFER_HOME'),'bin',...
  'mri_nu_correct.mni')), 'FREESURFER not accessible!')

%% 0-4. ANTs
if Job.UseAnts
  LD_LIBRARY_PATH = getenv('LD_LIBRARY_PATH');
  PATH_ = strsplit(LD_LIBRARY_PATH,':');
  ind = contains(upper(PATH_),'ANTS');
  assert(sum(ind), 'No ANTs library in LD_LIBRARY_PATH?!')
end

%% 0-5. SPM
assert(exist('spm', 'file'), 'SPM is not accessible!')
a = spm('Ver');
if ~strcmp(a(4:5),'12')
  error('Run this script on SPM12!');
end
spm('Defaults','fmri');
spm_jobman('initcfg');


%% ---T1w processing---
%% 1. Unified segmentation
% outputs:
% [1]  ${t1w}_seg8.mat : deformation field in native space & meta
% [2]  y_${t1w}.nii    : cosine functions in template space
% [3]  c?${t1w}.nii    : tissue prob maps in native space
% [4]  m${t1w}.nii     : bias-corrected t1w
% [5]  bm${t1w}.nii    : skull-stripped brain image
% [6]  wbm${t1w}.nii   : skull-stripped brain image registered in mni152

Job.FnameT1w = fullfile(Job.T1wDname, [Job.T1wPrefix,'.nii']);
Job.FnameT1wBm = fullfile(Job.T1wDname, ['bm',Job.T1wPrefix,'.nii']);
if not(isfile(Job.FnameT1wBm))
  tic; logthis('START: Unified segmentation of T1w\n')
  Job_ = struct('fname_t1w',Job.FnameT1w);
  Job_.norm = 1;
  switch (Job.Segment)
    case 'spm12'
      myspm_seg12(Job_, 'ss');
    case 'cat12'
      myspm_seg_cat12(Job_, 'ss');
    otherwise
      error('Job.Segment="%s" UNDEFINED METHOD', Job.Segment)
  end

  logthis('File created: '); ls(Job.FnameT1wBm)
  logthis('END: Unified segmentation of T1w: '); toc
end
Job.FnameT1wDeform = fullfile(Job.T1wDname, ['y_',Job.T1wPrefix,'.nii']);
assert(isfile(Job.FnameT1wDeform))


%% ANTS registration on skull-stripped T1w to skull-stripped MNI template
if Job.UseAnts
  % outputs:
  % [1]  ${t1w}_to_MNI152_T1_1mm_brain_stage2_Composite.h5: deformation field
  % [2]  ${t1w}_to_MNI152_T1_1mm_brain_stage2_fwdwarped.nii.gz: %{t1w} in MNI
  % [3]  ${t1w}_to_MNI152_T1_1mm_brain_stage2_invwarped.nii.gz: MNI in ${t1w}
  
  Job_ = [];
  Job_.fname_fixed = fullfile(getenv('FSLDIR'),'data','standard',...
    'MNI152_T1_1mm_brain.nii.gz');
  Job_.fname_moving = Job.FnameT1wBm;
  Job_.reg_stages = 2;
 
  Job.FnameT1wAntsDeform = fullfile(Job.T1wDname, ...
    ['bm',Job.T1wPrefix,'_to_MNI152_T1_1mm_brain_stage2_Composite.h5']);
  
  if not(isfile(Job.FnameT1wAntsDeform))
    logthis('START: ANTS: native T1w -> MNI152...\n'); tic
    myants_antsRegistration(Job_);
    logthis('END: ANTS: native T1w -> MNI152: '); toc
    logthis('File created: '); ls(Job.FnameT1wAntsDeform)
  end
  
end


%% ---EPI processing---

%% 2. Find meta-info, (crop, and) slice timing correction
% output:
% [1]  a${epi}.nii   : slice-timing corrected EPI
Job = myspm_tempproc(Job);


%% 3. unwarp+realign to MEAN IMAGE
% outputs:
% [1]  rp_a${epi}.txt   : six rigid-body motion parameters [mm & rad]
% [2]  a${epi}.mat      : [4x4xT] realign transform
% [3]  a${epi}_uw.mat   : unwarping meta data
% [4]  ua${epi}.nii     : unwarped/realigned image
% [5]  meanua${epi}.nii : mean image of [4]
Job_ = Job;
Job_.fname_epi = [Job.EpiDname, filesep, 'a',Job.EpiPrefix,'.nii'];
if exist('totalreadout_msec','var')
  Job_.totalreadout_msec = totalreadout_msec;
end
Job.FnameEpiMean = fullfile(Job.EpiDname,['meanua',Job.EpiPrefix,'.nii']);
if not(isfile(Job.FnameEpiMean))
  myspm_unwarp(Job_);
  logthis('File created: '); ls(Job.FnameEpiMean)
end

%% Intensity-bias correction of EPI (for stable coregistration)
Job.FnameEpiUnbiased = fullfile(Job.EpiDname, ['mmeanua',Job.EpiPrefix,'.nii']); % bias-corrected EPI
try
  if not(isfile(Job.FnameEpiUnbiased))
    system(['mri_nu_correct.mni --i ',Job.FnameEpiMean,' --o ',Job.FnameEpiUnbiased])
    logthis('File created: '); ls(Job.FnameEpiUnbiased)
  end
catch
  system(['ln -sf ',Job.FnameEpiMean,' ',Job.FnameEpiUnbiased])
  logthis('Link is created: '); ls(Job.FnameEpiUnbiased)
end


%% 4. COMPCOR regressors
Job_ = Job;
Job_.FnameEpi = fullfile(Job.EpiDname, ['ua',Job.EpiPrefix,'.nii']);
Job.FnameEigvec = fullfile(Job.EpiDname, ...
  sprintf('ua%s_n%i_b%.2f-%.2f_eigenvec.txt', ...
  Job.EpiPrefix, Job_.nPcs, Job_.Bpf));
if not(isfile(Job.FnameEigvec))
  myy_compcor(Job_);
  logthis('File created: '); ls(Job.FnameEigvec)
end


%% 5. Coregistration of EPI to native T1w
% outputs:
% - if SPM12 used:
% [1]  ua${epi}.nii       : header-modified EPI images
% [2]  rmmeanua${epi}.nii : resampled mean EPI image (for quality check)
%
% - if ANTs used:
% [1]  mmeanuaa${epi}_to_bm${t1w}_stage0_Composite.h5
%      : rigid-body transform (by ANTs)

if not(Job.UseAnts)
  % OUTPUT? resampled volume#1 for sanity check
  FnameOut = [Job.EpiDname, filesep, 'rmmeanua',Job.EpiPrefix,Ext];
  Job_ = Job;
  Job_.FnameEpi = [Job.EpiDname, filesep, 'ua',Job.EpiPrefix,Ext];
  if ~isfile(FnameOut)
    logthis('Coreg EPI to T1w using SPM12..\n')
    myspm_coreg_hdr(Job_);
    logthis('File created: '); ls(FnameOut)
  end
  
else  
  % RIGID from EPI to T1w
  Job_ = [];
  Job_.fname_moving = Job.FnameEpiUnbiased;  % bias-corrected EPI
  Job_.fname_fixed  = Job.FnameT1wBm;  % bias-corrected (skull-stripped) T1w
  Job_.reg_stages = 0;
  FnameOut = fullfile(Job.EpiDname, ...
    ['mmeanua',Job.EpiPrefix,'_to_bm',Job.T1wPrefix,'_stage0_Composite.h5']);
  if ~isfile(FnameOut)
    logthis('Coreg EPI to T1w using ANTs..\n')
    myants_antsRegistration(Job_);
    logthis('File created: '); ls(FnameOut)
  end
  
end



%% 6. Apply forward deformation on EPI for registration into MNI152
% outputs:
% - if SPM12 used:
% [1]  wua${epi}.nii    : resampled EPI images in MNI152 (by SPM12)
%
% - if ANTs used:
% [1]  ua${epi}_mni.nii : resampled EPI images in MNI152 (by ANTs)

if not(Job.UseAnts)
  %% EPI normalization using SPM12
  Job.EpiNorm = [Job.EpiDname, filesep, 'wua',Job.EpiPrefix,Ext];
  if ~isfile(Job.EpiNorm)
    logthis('Resampling EPI in MNI152 using SPM12..\n')
    myspm_norm(struct(...
      'fname_moving',[Job.EpiDname, filesep, 'ua',Job.EpiPrefix,Ext], ...
      'fname_deform',Job.FnameT1wDeform, 'vox_mm',Job.Vox_mm, 'interp',4, ...
      'bbox_mm',Job.Bbox_mm, 'sanitycheck',1));
    logthis('File created: '); ls(Job.EpiNorm)
  end

else
  %% EPI normalization using ANTs

  % CREATE FUNCTIONAL REFERENCE AT GIVEN VOXE RESOLUTION
  Job.FnameMniLow = fullfile(Job.EpiDname, 'mni_funcref.nii');
  if ~isfile(Job.FnameMniLow)
    % reslice
    fname_mni = fullfile(getenv('FSLDIR'), ...
      'data', 'standard', 'MNI152_T1_1mm_brain.nii.gz');
    % NOTE: both MNI152_T1_1mm and MNI152_T1_1mm_brain are in [182x218x182]
    system(['mri_convert -vs ',num2str(Job.Vox_mm),...
      ' ',fname_mni,' ',Job.FnameMniLow])
    % bounding-box
    if isfield(Job,'Bbox_mm') && ~isempty(Job.Bbox_mm)
      myspm_boundingbox(Job.FnameMniLow, Job.Bbox_mm, Job.FnameMniLow);
    else
      myspm_boundingbox(Job.FnameMniLow, 'canon', Job.FnameMniLow);
    end
    logthis('File created: '); ls(Job.FnameMniLow)
  end
  
  % COMBINE TRANSFORMS
  FnWarp = [Job.T1wDname, filesep, ...
    'mmeanua',Job.EpiPrefix,'_to_mni_Warping.nii.gz'];
  FnReg_epi_to_t1w = [Job.EpiDname, filesep, ...
    'mmeanua',Job.EpiPrefix,'_to_bm',Job.T1wPrefix,'_stage0_Composite.h5'];
  FnReg_t1w_to_mni = [Job.T1wDname, filesep, ...
    'bm',Job.T1wPrefix,'_to_MNI152_T1_1mm_brain_stage2_Composite.h5'];
  Job_ = struct('FnameOut',FnWarp, ...
    'FnameFixed',Job.FnameMniLow, ...
    'Transforms',{{FnReg_epi_to_t1w,0; FnReg_t1w_to_mni,0}});
  if ~isfile(FnWarp)
    logthis('Combining two level transforms in ANTs..\n')
    myants_combinetransforms(Job_);
    logthis('File created: '); ls(FnWarp)
  end
  
  % APPLY ON TIMESERIES
  Job.EpiNorm = [Job.EpiDname, filesep, 'xua',Job.EpiPrefix,Ext];
  Job_ = struct(...
    'FnameMoving',[Job.EpiDname, filesep, 'ua',Job.EpiPrefix,Ext], ...
    'FnameFixed',Job.FnameMniLow, 'Transforms',{{FnWarp,0}}, ...
    'FnameOut',Job.EpiNorm);
  if ~isfile(Job.EpiNorm)
    logthis('Resampling EPI in MNI152 using ANTs..\n')
    myants_antsApplyTransformsTimeseries(Job_);
    logthis('File created: '); ls(Job.EpiNorm)
  end
  
end


%% (7). Smoothing
if not(isempty(Job.Fwhm_mm))
  FnameOut = strrep(Job.EpiNorm, [Job.EpiDname, filesep], ...
    [Job.EpiDname, filesep, 's',num2str(Job.Fwhm_mm(1))]);
  if ~isfile(FnameOut)
    logthis('Spatial smoothing: FWHM=[%g] mm\n', Job.Fwhm_mm);
    myspm_smooth(struct('fname',Job.EpiNorm, 'fwhm_mm',Job.Fwhm_mm));
    logthis('File created: '); ls(FnameOut)
  end
end


%% Close the logfile
logfile(Job.FnameLog, Job, false, tStart);

end
