function JOB = myfsl_dwipreproc_woeddy (JOB)
%
% JOB requires:
% .fn_dwi    (raw data)
% (.fn_bval) if not "${fn_dwi}.bval"
% (.fn_bvec) if not "${fn_dwi}.bvec"
% .fn_t1w
%
% (cc) 2019, sgKIM  solleo@gmail.com

ls(JOB.fn_dwi)
[p1,f1,e1] = myfileparts(JOB.fn_dwi);
if ~strcmp(e1,'.nii')
  error('just give me .nii file! (it''s faster!)')
end
if ~isfield(JOB,'fn_bval')
  JOB.fn_bval = [p1,'/',f1,'.bval'];
end
ls(JOB.fn_bval)
if ~isfield(JOB,'fn_bvec')
  JOB.fn_bvec = [p1,'/',f1,'.bvec'];
end
ls(JOB.fn_bvec)

% Create nodiff.nii by averaging all images with bval=0
JOB.fn_nodiff = [p1,'/nodiff',e1];
nii = load_untouch_nii(JOB.fn_dwi);
bval = load(JOB.fn_bval);
nii.img = mean(nii.img(:,:,:,bval==0),4,'native');
nii.hdr.dime.dim(5) = 1;
nii.hdr.dime.pixdim(5) = 0;
save_untouch_nii(nii, JOB.fn_nodiff)
ls(JOB.fn_nodiff)

% Skull-stripping nodiff
myspm_seg12(JOB.fn_nodiff,'ss');

% Create brainmask
[p2,f2,e2] = myfileparts(JOB.fn_nodiff);
setenv('FSLOUTPUTTYPE','NIFTI');
fn_mask=[p2,'/bm',f2,'_mask',e2];
unix(['fslmaths ',p2,'/bm',f2,e2,' -bin -kernel 3D -ero -dilM ',fn_mask]);

% fit a tensor model:
unix(['dtifit ',...
  ' -k ',JOB.fn_dwi,' -o ',f1,...
  ' -m ',fn_mask,' -r ',f1,'.bvec -b ',f1,'.bval'])

% Denoising (smoothing along 4-th dimension)
addpath ~/Dropbox/MATLABfileExchange/DWIDenoising/
fname = JOB.fn_dwi;
method = 6; % local PCA
suffix = '_lpca';
lpca_prefix = [f1,'_lpca'];
JOB.fn_lpcaout = [lpca_prefix,e1];
dash = repmat('-',[1 72]);
if ~exist(JOB.fn_lpcaout, 'file')
  fprintf('%s\n[%s] Local-PCA (ETA=30 min)...\n',dash,datestr(now,31))
  %if ~strcmp(get(0,'diaryfile'),'diary'), diary('off');diary('on'); end
  tic
  MainDWIDenoising_cmd(fname, method, suffix)
  toc % Elapsed time is 546.644039 seconds. (9 min) on linux
end

% copy bval & bvec:
unix(['cp ',f1,'.bval ',lpca_prefix,'.bval']);
unix(['cp ',f1,'.bvec ',lpca_prefix,'.bvec']);

% fit a tensor model:
unix(['dtifit ',...
  ' -k ',JOB.fn_lpcaout,' -o ',lpca_prefix,...
  ' -m ',fn_mask,' -r ',lpca_prefix,'.bvec -b ',lpca_prefix,'.bval'])

% visualize improvement:
fa_orig = load_untouch_nii([f1,'_FA',e1]);
v1_orig = load_untouch_nii([f1,'_V1',e1]);
fa_lpca = load_untouch_nii([lpca_prefix,'_FA',e1]);
v1_lpca = load_untouch_nii([lpca_prefix,'_V1',e1]);
fn_gif = [lpca_prefix,'_vs_orig_v1.gif'];
figure('color','k', 'position',[642 439 1279 333],'visible','off')
Z = round(linspace(1, fa_orig.hdr.dime.dim(4), 19));
Z = Z(6:10);
Title = {JOB.fn_dwi, JOB.fn_lpcaout};
ax = axeslayout([1 5],[0 0 0.05 0],[0 0 0 0]);
I = {v1_orig.img, v1_lpca.img};
J = {fa_orig.img, fa_lpca.img};
delay_sec = [2 2];
for i=1:2
  clf
  for j=1:5
    axespos(ax, j)
    rgbslice = permute(abs(squeeze(I{i}(:,:,Z(j),:))),[2 1 3]);
    image(rgbslice .* 1.5 .* zeroone(repmat(J{i}(:,:,Z(j))',[1 1 3])))
    axis off; axis image; set(gca,'ydir','nor')
    axis([
      find(sum(sum(rgbslice,3),1),1,'first')
      find(sum(sum(rgbslice,3),1),1,'last')
      find(sum(sum(rgbslice,3),2),1,'first')
      find(sum(sum(rgbslice,3),2),1,'last')]);
    if j==3
      title(Title{i},'color','w', 'interp','none','fontweight','normal');
    end
  end
  gifani(gcf, fn_gif, i, delay_sec(i));
end

% GZIP:
if exist([lpca_prefix,'.nii'],'file') && ~exist([lpca_prefix,'.nii.gz'],'file')
  unix(['gzip ',lpca_prefix,'.nii'])
elseif exist([lpca_prefix,'.nii'],'file') && exist([lpca_prefix,'.nii.gz'],'file')
  delete([lpca_prefix,'.nii'])
end

end
