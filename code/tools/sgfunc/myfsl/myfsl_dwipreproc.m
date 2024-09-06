function JOB = myfsl_dwipreproc (JOB)
%
% JOB requires:
% .fn_dwi    (raw data)
% (.fn_bval) if not "${fn_dwi}.bval"
% (.fn_bvec) if not "${fn_dwi}.bvec"
% .fn_pe_pair
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

% TOPUP: Unwarp nodiff
dash = repmat('-',[1 60]);
fprintf('%s\n[%s] Topup (ETA=10 min)...\n',dash,datestr(now,31))
%if ~strcmp(get(0,'diaryfile'),'diary'), diary('off');diary('on'); end
JOB2 = JOB;
JOB2.fn_to_unwarp = JOB.fn_nodiff;
JOB2 = myfsl_topup_siemens(JOB2);

% Skull-stripping nodiff
myspm_seg12(JOB2.fn_unwarped,'ss');

% Create brainmask
[p2,f2,e2] = myfileparts(JOB2.fn_unwarped);
setenv('FSLOUTPUTTYPE','NIFTI');
fn_mask=[p2,'/bm',f2,'_mask',e2];
unix(['fslmaths ',p2,'/bm',f2,e2,' -bin -kernel 3D -ero -dilM ',fn_mask]);

% EDDY: eddy current correction + susceptability distortion unwarping + 
% signal drop out due to subject motion correction + outlier slice correction
fprintf('%s\n[%s] Eddy (ETA=3 hr)...\n',dash,datestr(now,31))
%if ~strcmp(get(0,'diaryfile'),'diary'), diary('off');diary('on'); end
fn_index = [p1,'/',f1,'_index.txt'];
fid = fopen(fn_index,'w');
for i=1:numel(bval); fprintf(fid,'1\n'); end; fclose(fid);
if ~exist([getenv('FSLDIR'),'/bin/eddy_openmp'],'file')
  eddy_cmd = 'eddy'; % this is depricated in FSL 5.0.11
else
  setenv('OMP_NUM_THREADS','16')
  eddy_cmd = 'eddy_openmp';
end
eddy_prefix = [p1,'/',f1,'_eddy'];
JOB.fn_eddyout = [eddy_prefix,e1];
if ~exist(JOB.fn_eddyout,'file')
  tic
  myunix([eddy_cmd, ...
    ' --imain=',JOB.fn_dwi,' --mask=',fn_mask, ...
    ' --bvecs=',JOB.fn_bvec,' --bvals=',JOB.fn_bval, ...
    ' --index=',fn_index,' --acqp=',JOB2.fn_acq, ...
    ' --out=',eddy_prefix,' --topup=',JOB2.prefix_out,' --repol']);
  toc 
  % Elapsed time is 8491.465967 seconds (2.3 hr on my macbook)
  % Elapsed time is 1590.460446 seconds (26 min on O-lab linux)
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

% Denoising (smoothing along 4-th dimension)
addpath ~/Dropbox/MATLABfileExchange/DWIDenoising/
fname = JOB.fn_eddyout;
method = 6; % local PCA
suffix = '_lpca';
lpca_prefix = [eddy_prefix,'_lpca'];
JOB.fn_lpcaout = [lpca_prefix,e1];
if ~exist(JOB.fn_lpcaout, 'file')
    fprintf('%s\n[%s] Local-PCA (ETA=30 min)...\n',dash,datestr(now,31))
    %if ~strcmp(get(0,'diaryfile'),'diary'), diary('off');diary('on'); end
    tic
    MainDWIDenoising_cmd(fname, method, suffix)
    toc % Elapsed time is 546.644039 seconds. (9 min) on linux
    
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
    figure('color','k', 'position',[642 439 1279 333],'visible','off')
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
end

% Recreate nodiff from eddy+lpca because this image can be used for
% accurate transforms between modalities.

nii = load_untouch_nii('dwi_eddy_lpca.nii');
bval = load('dwi_eddy_lpca.bval');
nii.img = mean(nii.img(:,:,:,bval==0),4);
nii.hdr.dime.dim(1) = 3;
nii.hdr.dime.dim(5) = 1;
save_untouch_nii(nii, 'nodiff_eddy_lpca.nii');
myspm_seg12('nodiff_eddy_lpca.nii','ss')

end
