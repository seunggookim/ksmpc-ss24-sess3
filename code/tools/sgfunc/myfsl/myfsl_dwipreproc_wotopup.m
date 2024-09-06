function JOB = myfsl_dwipreproc_wotopup (JOB)
%
% JOB requires:
% .fn_dwi    (raw data)
% (.fn_bval) if not "${fn_dwi}.bval"
% (.fn_bvec) if not "${fn_dwi}.bvec"
% .fn_t1w
%
% Ref: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FDT/UserGuide
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
if ~isfield(JOB,'fn_bvec')
  JOB.fn_bvec = [p1,'/',f1,'.bvec'];
end

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
myspm_seg12(JOB.fn_nodiff,'ss')

% Create brainmask
[p2,f2,e2] = myfileparts(JOB.fn_nodiff);
setenv('FSLOUTPUTTYPE','NIFTI');
fn_mask=[p2,'/bm',f2,'_mask',e2];
unix(['fslmaths ',p2,'/bm',f2,e2,' -bin -kernel 3D -ero -dilM ',fn_mask]);

% Eddy current correction + susceptability distortion unwarping + signal
% drop out due to subject motion correction + outlier slice correction
fn_index = [p1,'/',f1,'_index.txt'];
fid = fopen(fn_index,'w');
for i=1:numel(bval); fprintf(fid,'1\n'); end; fclose(fid);


%% READ HEADERS for acquisition parameters
fn_acq = [p1,'/',f1,'_acq.txt'];
JOB.fn_acq = fn_acq;
if ~exist(fn_acq, 'file')
mtx = {};
n=1;
for i=1
  json{i} = jsondecode(fileread([p1,'/',f1,'.json']));
  switch json{i}.PhaseEncodingDirection(1)
    case 'i' % PhaseEncodingAxis: 'i' or 'i-' : signs are only available for Siemens scanners
      mtx{i} = [1 0 0];
    case 'j' 
      mtx{i} = [0 1 0];
    case 'k'
      mtx{i} = [0 0 1];
  end
  if numel(json{i}.PhaseEncodingDirection)==2 ...
      && strcmp(json{i}.PhaseEncodingDirection(2),'-')
    mtx{i} = -mtx{i};
  end
end
M=[
  repmat([mtx{1} json{1}.TotalReadoutTime],[n 1])
  ];
dlmwrite(fn_acq,M,'\t')
end
%%

[~,stderr]=unix('which eddy_openmp');
if isempty(stderr)
  eddy_cmd = 'eddy'; % this is depricated in FSL 5.0.11
else
  setenv('OMP_NUM_THREADS','16')
  eddy_cmd = 'eddy_openmp';
end
eddy_prefix = [p1,'/',f1,'_eddy'];
JOB.fn_eddyout = [eddy_prefix,e1];
if ~exist(JOB.fn_eddyout,'file')
  tic
  unix([eddy_cmd, ...
    ' --imain=',JOB.fn_dwi,' --mask=',fn_mask, ...
    ' --bvecs=',JOB.fn_bvec,' --bvals=',JOB.fn_bval, ...
    ' --index=',fn_index,' --acqp=',JOB.fn_acq, ...
    ' --out=',f1,'_eddy'])
  toc % Elapsed time is 8491.465967 seconds (2.3 hr on my macbook)
end
% Rename rotated bvecs and copy bvals for the eddy output
unix(['cp ',JOB.fn_bval,' ',eddy_prefix,'.bval']);
unix(['cp ',p1,'/',f1,'_eddy.eddy_rotated_bvecs ',eddy_prefix,'.bvec']);

% fit a tensor model:
unix(['dtifit ',...
  ' -k ',JOB.fn_eddyout,' -o ',eddy_prefix,...
  ' -m ',fn_mask,' -r ',eddy_prefix,'.bvec -b ',eddy_prefix,'.bval'])

% Visualize outlier corrections:
fn_tmp = tempname;
unix(['sed "1d;$d" ',p1,'/',f1,'_eddy.eddy_outlier_map > ',fn_tmp]);
outmap = load(fn_tmp);
[iScans,iSlices] = find(outmap);
if ~isempty(iScans)
  figure('position',[1190 776 731 270],'visible','off')
  ax = axeslayout([1 3],'tight','tight');
  for j=1:numel(iScans)
    rawVol  = load_untouch_nii(JOB.fn_dwi, iScans(j));
    eddyVol = load_untouch_nii(JOB.fn_eddyout, iScans(j));
    Im = {rawVol.img(:,:,iSlices(j))',eddyVol.img(:,:,iSlices(j))'};
    clf
    axespos(ax,1);
    imagesc(Im{1});
    set(gca,'ydir','nor'); axis image;
    title(sprintf('Raw: scan#=%i, slice#=%i',iScans(j),iSlices(j)))
    
    axespos(ax,2);
    imagesc(Im{2});
    set(gca,'ydir','nor'); axis image;
    title('Eddy corrected')
    xlabel(JOB.fn_eddyout,'interp','none')
    
    axespos(ax,3);
    imagesc(Im{1}-Im{2});
    set(gca,'ydir','nor'); axis image;
    title('Difference')
    
    colormap gray
    export_fig([p1,'/',f1,'_eddy.eddy_outlier',num2str(j),'.png'],'-r150')
  end
  close(gcf)
end

% delete "eddy_outlier_free" file (This is corrected ONLY for outliers but
% NOT for all the other artifacts. the final output is "outputprefix.nii")
% just confusing for me)
unix(['rm -f ',p1,'/',f1,'_eddy.eddy_outlier_free_data',e1])

% (recreate nodiff? it's possible that bval=0 images could be also corrected
% but it's only used for brainmask and not for the diffusion modeling. so
% not really?)

% Denoising (smoothing along 4-th dimension)
addpath ~/Dropbox/MATLABfileExchange/DWIDenoising/
fname = JOB.fn_eddyout;
method = 6; % local PCA
suffix = '_lpca';
lpca_prefix = [eddy_prefix,'_lpca'];
JOB.fn_lpcaout = [lpca_prefix,e1];
if ~exist(JOB.fn_lpcaout, 'file')
  try
    tic
    MainDWIDenoising_cmd(fname, method, suffix)
    toc % Elapsed time is 1856.644039 seconds.
    
    % copy bval & bvec:
    unix(['cp ',eddy_prefix,'.bval ',lpca_prefix,'.bval']);
    unix(['cp ',eddy_prefix,'.bvec ',lpca_prefix,'.bvec']);
    
    % fit a tensor model:
    unix(['dtifit ',...
      ' -k ',JOB.fn_lpcaout,' -o ',lpca_prefix,...
      ' -m ',fn_mask,' -r ',lpca_prefix,'.bvec -b ',lpca_prefix,'.bval'])
    
    % visualize improvement:
    fa_eddy = load_untouch_nii([eddy_prefix,'_FA',e1]);
    v1_eddy = load_untouch_nii([eddy_prefix,'_V1',e1]);
    fa_lpca = load_untouch_nii([lpca_prefix,'_FA',e1]);
    v1_lpca = load_untouch_nii([lpca_prefix,'_V1',e1]);
    fn_gif = [lpca_prefix,'_vs_eddy_v1.gif'];
    figure('color','k', 'position',[642 439 1279 333],'visible','on')
    Z = round(linspace(1, fa_eddy.hdr.dime.dim(4), 19));
    Z = Z(6:10);
    Title = {JOB.fn_eddyout, JOB.fn_lpcaout};
    ax = axeslayout([1 5],[0 0 0.05 0],[0 0 0 0]);
    I = {v1_eddy.img, v1_lpca.img};
    J = {fa_eddy.img, fa_lpca.img};
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
  catch
    warning('Could not find DWIDenoising path: denoising not applied')
  end
end

end