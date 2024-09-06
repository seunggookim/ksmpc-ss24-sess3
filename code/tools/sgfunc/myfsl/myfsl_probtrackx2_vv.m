function myfsl_probtrackx2_vv(JOB)
% directory convention = "imgDB/20548/fdt/data"

if ~exist('JOB','var'), JOB=[]; end
if ~isfield(JOB,'wm_thr'), JOB.wm_thr=0.99; end
if ~isfield(JOB,'nsamples'), JOB.nsamples=5000; end
if ~isfield(JOB,'dir_fdt'), JOB.dir_fdt=pwd; end

cd(JOB.dir_fdt)
if ~exist('nodif.nii','file')
 bvals=load('bvals');
 nii=load_untouch_nii('data.nii.gz');
 nii2=nii;
 nii2.img=int16(mean(nii.img(:,:,:,bvals==0),4));
 nii2.hdr.dime.dim(1)=3;
 nii2.hdr.dime.dim(5)=1;
 save_untouch_nii(nii2,'nodif.nii');
end
%% Create brain mask, bias-corrected EPI (for coreg)
if ~exist('nodif_brain_mask.nii.gz','file')
 JOB=[];
 JOB.native=[1 1 1 0 0 0];
 JOB.mw=zeros(1,6);
 JOB.fname_t1w='nodif.nii';
 myspm_seg12(JOB);
 setenv('FSLOUTPUTTYPE','NIFTI_GZ')
 unix('fslmaths bmnodif.nii -bin nodif_brain_mask')
end
%% BEDPOSTX
if ~exist('merged_f1samples.nii.gz','file')
 unix('bedpostx .')
end
%% Bring GM, WM masks based on T1w (Hmmph!)
if ~exist('gm_mask.nii.gz','file') || ~exist('wm_mask.nii.gz','file')
 unix(['cp ../../spm12/bmt1w.nii .']);
 unix(['cp ../../spm12/c1t1w.nii .']);
 unix(['cp ../../spm12/c2t1w.nii .']);
 JOB=[];
 JOB.interp=1;
 JOB.prefix='d';
 JOB.fname_moving='bmt1w.nii';
 JOB.fname_fixed='bmnodif.nii';
 JOB.fname_others={'c1t1w.nii','c2t1w.nii'};
 myspm_coreg(JOB);
 unix('rm bmt1w.nii');
 unix('rm c1t1w.nii');
 unix('rm c2t1w.nii');
 setenv('FSLOUTPUTTYPE','NIFTI_GZ')
 unix('fslmaths dc1t1w.nii -thr .80 -bin gm_mask');
 unix(['fslmaths dc2t1w.nii -thr ',num2str(JOB.wm_thr),' -bin wm_mask']);
end

%% Matrix3: Seeding from all white matter voxels,
% bidirectionally propogate, until reach targets (all GM voxels?)
% constraining tracking by cortex surfaces (to prevent a spurious
% tracking that connects adjacent gyri)

dir_fdt=[JOB.dir_fdt];
dir_bpx=[dir_fdt,'.bedpostX'];
dir_ptx=[dir_fdt,'/net.wm',num2str(JOB.wm_thr*100),'.n',num2str(JOB.nsamples)];

if ~exist([dir_fdt,'/net.wb/fdt_matrix3.dot'],'file')
 unix(['probtrackx2 ', ...
  ' --samples=',dir_bpx,'/merged ', ...
  ' --mask=',dir_fdt,'/nodif_brain_mask.nii.gz ', ...
  ' --seed=',dir_fdt,'/wm_mask.nii.gz ', ...
  ' --omatrix3 ', ...
  ' --target3=',dir_fdt,'/gm_mask.nii.gz ', ...
  ' --loopcheck --opd --verbose=1', ...
  ' --nsamples=',num2str(JOB.nsamples), ...
  ' --distthresh=3 ', ...
  ' --dir=',dir_ptx, ...
  ' --forcedir', ...
  ' --out=all '])
end

end