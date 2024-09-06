function myspm_postmeica(JOB)
% JOB requires:
% .fname_epi
% .fname_t1w
% .TR_sec
% .vox_mm

if ~isfield(JOB,'overwrite'), JOB.overwrite=0; end
spm('Defaults','fmri');
spm_jobman('initcfg');
%% ---T1w processing---
%% 1. Unified segmentation
% outputs:
% [1]  ${t1w}_seg8.mat : deformation field in native space & meta
% [2]  y_${t1w}.nii    : cosine functions in template space (wtf?)
% [3]  c?${t1w}.nii    : tissue prob maps in native space
% [4]  mwc?${t1w}.nii  : JD modulated TPM in template space
% [5]  m${t1w}.nii     : bias-corrected t1w
% [6]  bm${t1w}.nii    : skull-stripped brain image
[p2,f2,e2]=myfileparts(JOB.fname_t1w);
JOB.fname_t1w=[p2,'/',f2,e2];
dir_tpm=[spm('dir'),'/tpm'];
matlabbatch={};
preproc=[];
preproc.channel.vols = {[JOB.fname_t1w,',1']};
preproc.channel.biasreg = 0.001;
preproc.channel.biasfwhm = 60;
preproc.channel.write = [0 1];
preproc.tissue(1).tpm = {[dir_tpm,'/TPM.nii,1']};
preproc.tissue(1).ngaus = 1;
preproc.tissue(1).native = [1 1]; % [native dartel]
preproc.tissue(1).warped = [0 1];
preproc.tissue(2).tpm = {[dir_tpm,'/TPM.nii,2']};
preproc.tissue(2).ngaus = 1;
preproc.tissue(2).native = [1 1]; % [native dartel]
preproc.tissue(2).warped = [0 1];
preproc.tissue(3).tpm = {[dir_tpm,'/TPM.nii,3']};
preproc.tissue(3).ngaus = 2;
preproc.tissue(3).native = [1 0];
preproc.tissue(3).warped = [0 1];
preproc.tissue(4).tpm = {[dir_tpm,'/TPM.nii,4']};
preproc.tissue(4).ngaus = 3;
preproc.tissue(4).native = [0 0];
preproc.tissue(4).warped = [0 0];
preproc.tissue(5).tpm = {[dir_tpm,'/TPM.nii,5']};
preproc.tissue(5).ngaus = 4;
preproc.tissue(5).native = [0 0];
preproc.tissue(5).warped = [0 0];
preproc.tissue(6).tpm = {[dir_tpm,'/TPM.nii,6']};
preproc.tissue(6).ngaus = 2;
preproc.tissue(6).native = [0 0];
preproc.tissue(6).warped = [0 0];
preproc.warp.mrf = 1;
preproc.warp.cleanup = 1;
preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
preproc.warp.affreg = 'mni';
if isfield(JOB,'isestern') && JOB.isestern
  preproc.warp.affreg = 'eastern';
end
preproc.warp.fwhm = 0;
preproc.warp.samp = 3;
preproc.warp.write = [0 1];
matlabbatch{1}.spm.spatial.preproc = preproc;
fname_in = JOB.fname_t1w;
ls(fname_in);
fname_out = [p2,'/bm',f2,e2];
if ~exist(fname_out,'file') || JOB.overwrite
  disp('[1] Unified segmentation of T1w..');
  fname_matlabbatch=[p2,'/spm12_fmriprep1.mat'];
  save(fname_matlabbatch,'matlabbatch');
  spm_jobman('run', matlabbatch);
  nii=cell(1,3);
  for c=1:3
    nii{c} = load_untouch_nii([p2,'/c',num2str(c),f2,e2]);
  end
  manat = load_untouch_nii([p2,'/m',f2,e2]);
  manat.img(~(nii{1}.img | nii{2}.img | nii{3}.img))=0;
  save_untouch_nii(manat, fname_out);
  clear manat nii
  isdone(fname_out,1);
end
%% 2. register the native brain to MNI152
% output:
% [1]  wb${t1w}.nii : transformed skull-stripped brain in template space
fname_in = [p2,'/bm',f2,e2];
ls(fname_in)
fname_def = [p2,'/y_',f2,e2];
ls(fname_def)
fname_out = [p2,'/wbm',f2,e2];
matlabbatch={};
write1=[];
write1.subj.resample = {fname_in};
write1.subj.def = {fname_def};
write1.woptions.bb = [-78 -112 -70; 78 76 85];
write1.woptions.vox = [1 1 1];
write1.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write = write1;
if ~exist(fname_out,'file') || JOB.overwrite
  disp('[2] Registering T1w to MNI..');
  fname_matlabbatch=[p2,'/spm12_fmriprep2.mat'];
  save(fname_matlabbatch,'matlabbatch');
  spm_jobman('run', matlabbatch);
  myspm_denan(fname_out,1);
  isdone(fname_out,2);
end
%% ---EPI processing---
% Find filenames
if ~iscell(JOB.fname_epi)
  JOB.fname_epi = {JOB.fname_epi};
end
for j = 1:numel(JOB.fname_epi)
  fname_epi = JOB.fname_epi{j};
  [p1,f1,e1] = myfileparts(fname_epi);
  fname_epi = [p1,'/',f1,e1];
  text1=['[Session ',num2str(j),']:'];
  disp([text1,repmat('-',1,60-numel(text1))]);
  ls(fname_epi)
  hdr = load_nii_hdr(fname_epi);
  NumFrames = hdr.dime.dim(5);
  disp(['> Number of frames = ',num2str(NumFrames)]);
  if ~isfield(JOB,'TR_sec')
    JOB.TR_sec=hdr.dime.pixdim(5);
    if ~JOB.TR_sec
      error(['TR=0 in the file header; Enter JOB.TR_sec !']);
    end
  end
  disp(['> TR = ',num2str(JOB.TR_sec),' sec']);
  %% 1. Compute mean
  setenv('FSLOUTPUTTYPE','NIFTI')
  unix(['fslmaths ',fname_epi,' -Tmean ',p1,'/mean',f1,e1]);
  %% 2. Intensity-bias correction of EPI (for better coregistration)
  fname_epi_unbiased = [p1,'/mmean',f1,e1]; % bias-corrected EPI
  if ~exist(fname_epi_unbiased,'file')
    unix(['mri_nu_correct.mni --i ',p1,'/mean',f1,e1,' --o ',fname_epi_unbiased]);
  end
  %% 3. Coregistration of EPI to native T1w
  % outputs: <modifying transform matrices in headers>
  % [1]  mmean${epi}.nii
  % [2]  ${epi}.nii        : "other" images, header modified
  % [3]  ${epi}.mat        : "other" images, header modified
  estimate1=[];
  fname_t1w_brain = [p2,'/bm',f2,e2]; % bias-corrected T1w
  estimate1.ref{1}    = fname_t1w_brain;
  ls(fname_t1w_brain)
  fname_epi_unbiased = [p1,'/mmean',f1,e1]; % bias-corrected EPI
  estimate1.source{1} = fname_epi_unbiased ;
  ls(fname_epi_unbiased)
  % ------------------here you specified OTHER imgaes-----------------
  fname_epi_ua = [p1,'/',f1,e1];
  ls(fname_epi_ua);
  estimate1.other{1} = fname_epi_ua; % NIFTI-format have one xfm for all volumes
  % ------------------------------------------------------------------
  estimate1.eoptions.cost_fun = 'nmi';
  estimate1.eoptions.sep = [4 2];
  estimate1.eoptions.tol = ...
    [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
  estimate1.eoptions.fwhm = [7 7];
  matlabbatch={};
  matlabbatch{1}.spm.spatial.coreg.estimate = estimate1;
  fname_out=[p1,'/',f1,'.mat'];
  if ~exist(fname_out,'file') || JOB.overwrite
    disp('[3] Coregistration of EPI to native T1w..');
    fname_matlabbatch=[p1,'/spm12_fmriprep6.mat'];
    save(fname_matlabbatch,'matlabbatch');
    spm_jobman('run', matlabbatch);
    isdone(fname_out,3);
  end
  %% 4. resample the first volume (to check coregistration quality)
  [p4,f4,e4]=myfileparts(fname_epi_unbiased);
  if ~exist([p4,'/r',f4,e4],'file')  || JOB.overwrite
    matlabbatch={};
    matlabbatch{1}.spm.spatial.coreg.write.ref{1}    = fname_t1w_brain;
    matlabbatch{1}.spm.spatial.coreg.write.source{1} = [fname_epi_unbiased,',1'];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
    spm_jobman('run', matlabbatch);
    fname_epi_in_t1w=[p1,'/rmmean',f1,e1];
    myspm_denan(fname_epi_in_t1w,1);
    unix(['slices ',fname_epi_in_t1w,' ',fname_t1w_brain, ...
      ' -o ',p1,'/rmmean',f1,'_in_',f2,'.gif'])
  end
  %% 5. Apply forward deformation on EPI for registration into MNI152
  matlabbatch={};
  write1=[];
  write1.subj.def{1} = fname_def;
  write1.subj.resample{1} = fname_epi_ua;
  write1.woptions.bb = [-78 -112 -70; 78 76 85];
  if ~isfield(JOB,'vox_mm')
    JOB.vox_mm = hdr.dime.pixdim(2:4);
    disp(['Resampling at [',num2str(JOB.vox_mm),'] mm in MNI152 space..'])
  end
  if numel(JOB.vox_mm) == 1
    JOB.vox_mm=[1 1 1]*JOB.vox_mm;
  end
  write1.woptions.vox = JOB.vox_mm;
  write1.woptions.interp = 4;
  matlabbatch{1}.spm.spatial.normalise.write = write1;
  fname_out=[p1,'/w',f1,e1];
  if ~exist(fname_out,'file') || JOB.overwrite
    disp('[5] Registration of EPI to MNI152..');
    fname_matlabbatch=[p1,'/spm12_fmriprep8.mat'];
    save(fname_matlabbatch,'matlabbatch');
    spm_jobman('run', matlabbatch);
    %myspm_denan(fname_out,1);
    addTR(fname_out,JOB.TR_sec);
    isdone(fname_out,5);
    
    % visual inspection
    fname_gif=[p1,'/w',f1,'_in_mni152.gif'];
    fname_epi_in_mni1=[p1,'/w',f1,'1.nii.gz'];
    setenv('FSLOUTPUTTYPE','NIFTI_GZ');
    unix(['fslroi ',p1,'/w',f1,e1,' ',fname_epi_in_mni1,' 0 1']);
    myspm_denan(fname_epi_in_mni1);
    ls(fname_epi_in_mni1);
    fname_mni=[getenv('FSLDIR'),'/data/standard/MNI152_T1_2mm.nii.gz'];
    fname_mni1=[p1,'/mni152.nii.gz'];
    unix(['mri_convert --like ',fname_epi_in_mni1,' ',fname_mni,' ',fname_mni1]);
    ls(fname_mni1)
    unix(['slices ',fname_epi_in_mni1,' ',fname_mni1,' -o ',fname_gif]);
    ls(fname_gif)
  end
  
  %% 6. Mask non-brain tissue voxels
  fname_out=[p1,'/bw',f1,e1];
  if ~exist(fname_out,'file')
    fname_epi_in_mni1=[p1,'/w',f1,'1.nii.gz'];
    fname_t1w_mni=[p2,'/wbm',f2,e2];
    fname_t1w_mni_low=[p2,'/wbm',f2,'_low',e2];
    fname_t1w_mni_low_mask=[p2,'/wbm',f2,'_low_mask',e2];
    unix(['mri_convert --like ',fname_epi_in_mni1,...
      ' ',fname_t1w_mni,' ',fname_t1w_mni_low]);
    setenv('FSLOUTPUTTYPE','NIFTI');
    unix(['fslmaths ',fname_t1w_mni_low,' -thr 10 -bin -s 3 -dilM -ero ',...
      '-thr 0.5 -bin ',fname_t1w_mni_low_mask]);
    unix(['fslmaths ',p1,'/w',f1,e1,' -mas ',fname_t1w_mni_low_mask,...
      ' ',p1,'/bw',f1,e1]);
    addTR([p1,'/bw',f1,e1],JOB.TR_sec);
    isdone(fname_out,6)
    unix(['slices ',p1,'/w',f1,e1,' ',fname_t1w_mni_low_mask,' -o ',...
      p1,'/w',f1,'_masked.gif']);
  end
end
end

%% Subroutines:
function isdone(fname_out,proci)
if exist(fname_out,'file')
  disp(['[',num2str(proci),'] done: ',fname_out])
else
  error(['[',num2str(proci),'] failed: ',fname_out,' was not created'])
end
end

function addTR(fname,TR)
nii=load_untouch_nii(fname);
nii.hdr.dime.pixdim(5)=TR;
save_untouch_nii(nii,fname);
end