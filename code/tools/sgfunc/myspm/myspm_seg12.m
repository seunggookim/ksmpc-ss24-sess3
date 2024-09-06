function JOB=myspm_seg12(JOB,runmode)
% JOB=myspm_seg12(JOB, [runmode]) runs SPM12 segmentation
%
% JOB requires:
%  .fname_t1w  [string]
% (.ismp2rage) [1x1] default=0 (if 1 use uncorrected input)
% (.iseastern) [1x1] default=0 (if 1 use East-Asian template for affine transform)
% (.mw)
% (.naitve)
% (.norm)
% (.runmode)   [string] 'SS'(mw=zeros(1,6), norm=1) | 
%
% if isMp2rage=1, then using __original image__ instead of intensity-bias
% corrected one (i.e., m*)
%
% (cc) 2017, sgKIM, solleo@gmail.com

ver = spm('version');
if ~strcmp(ver(4:5),'12'), error('This function is for SPM12!'); end
if nargin == 0, help(mfilename); return; end
if ~isstruct(JOB) && ischar(JOB), JOB=struct('fname_t1w',JOB); end
if ~isfield(JOB,'ismp2rage'), JOB.ismp2rage=0; end
if ~isfield(JOB,'iseastern'), JOB.iseastern=0; end
[p2,f2,e2]=myfileparts(JOB.fname_t1w);
dir_tpm=[spm('dir'),filesep,'tpm'];
if ~isfield(JOB,'norm'), JOB.norm=0; end
if exist('runmode','var')
  switch runmode
    case {'ss'}
      JOB.mw=zeros(1,6);
      JOB.norm=1;
  end
end

preproc=[];
preproc.channel.vols={[JOB.fname_t1w,',1']};
preproc.channel.biasreg=0.001;
preproc.channel.biasfwhm=60;
preproc.channel.write=[0 1];
preproc.tissue(1).tpm={[dir_tpm,filesep,'TPM.nii,1']};
preproc.tissue(1).ngaus=1;
preproc.tissue(1).native=[1 0];
preproc.tissue(1).warped=[0 1];
preproc.tissue(2).tpm={[dir_tpm,filesep,'TPM.nii,2']};
preproc.tissue(2).ngaus=1;
preproc.tissue(2).native=[1 0];
preproc.tissue(2).warped=[0 1];
preproc.tissue(3).tpm={[dir_tpm,filesep,'TPM.nii,3']};
preproc.tissue(3).ngaus=2;
preproc.tissue(3).native=[1 0];
preproc.tissue(3).warped=[0 1];
preproc.tissue(4).tpm={[dir_tpm,filesep,'TPM.nii,4']};
preproc.tissue(4).ngaus=3;
preproc.tissue(4).native=[0 0];
preproc.tissue(4).warped=[0 0];
preproc.tissue(5).tpm={[dir_tpm,filesep,'TPM.nii,5']};
preproc.tissue(5).ngaus=4;
preproc.tissue(5).native=[0 0];
preproc.tissue(5).warped=[0 0];
preproc.tissue(6).tpm={[dir_tpm,filesep,'TPM.nii,6']};
preproc.tissue(6).ngaus=2;
preproc.tissue(6).native=[0 0];
preproc.tissue(6).warped=[0 0];
if isfield(JOB,'mw')
  for c=1:6
    preproc.tissue(c).warped(2)=JOB.mw(c);
  end
end
if isfield(JOB,'native')
  for c=1:6
    preproc.tissue(c).native(1)=JOB.native(c);
  end
end
preproc.warp.mrf=1;
preproc.warp.cleanup=1;
preproc.warp.reg=[0 0.001 0.5 0.05 0.2];
if JOB.iseastern
  preproc.warp.affreg='eastern';
else
  preproc.warp.affreg='mni';
end
preproc.warp.fwhm=0;
preproc.warp.samp=3;
preproc.warp.write=[0 1];

matlabbatch={};
matlabbatch{1}.spm.spatial.preproc=preproc;
ls(JOB.fname_t1w);
if JOB.ismp2rage
  fname_out=[p2,filesep,'b',f2,e2];
else
  fname_out=[p2,filesep,'bm',f2,e2];
end
if ~exist(fname_out,'file')
  fname_matlabbatch=[p2,filesep,'spm12_seg.mat'];
  save(fname_matlabbatch,'matlabbatch');
  spm_jobman('initcfg');
  spm_jobman('run',matlabbatch);
  Y=cell(1,3);
  for c=1:3
    V=spm_vol_nifti([p2,filesep,'c',num2str(c),f2,e2]);
    [Y{c},~]=spm_read_vols(V);
  end
  if JOB.ismp2rage
    V=spm_vol_nifti([p2,filesep,'',f2,e2]);
  else
    V=spm_vol_nifti([p2,filesep,'m',f2,e2]);
  end
  manat_uncut=spm_read_vols(V);
  manat = manat_uncut;
  manat(~(Y{1}|Y{2}|Y{3}))=0;
  V.fname=fname_out;
  V.descrip='Skull stripped';
  spm_write_vol(V,manat);

  try
  % quality check:
  fname_png=[p2,filesep,'bm',f2,'.png'];
  slices(manat_uncut, struct('fname_png',fname_png), double(Y{1}|Y{2}|Y{3}))
  end
end

if JOB.norm
  [p3,f3,e3] = myfileparts(fname_out);
  if ~exist([p3,filesep,'w',f3,e3],'file')
    job1 = [];
    job1.fname_deform = [p2,filesep,'y_',f2,e2];
    job1.vox_mm = [1 1 1];
    job1.fname_moving = fname_out;
    myspm_norm(job1)
    fname_out=[p3,filesep,'w',f3,e3];
    ls(fname_out);
    
    try
    % quality check:
    fname_png=[p3,filesep,'w',f3,'_in_mni152.png'];
    fname_mni=[getenv('FSLDIR'),filesep,'data/standard/MNI152_T1_1mm_brain.nii.gz'];
    slices(fname_out, struct('fname_png',fname_png), fname_mni)
    end
  end
end

disp('> Done:');
ls(fname_out)
end
