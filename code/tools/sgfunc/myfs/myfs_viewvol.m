function [H, cfg] = myfs_viewvol( base, data, cfg )
%MYFS_VIEWVOL plots slices from given volumes
% [H, cfg] = myfs_viewvol( base, data, cfg )

% NOTE
%{
nifti.vox2ras is 0-based!
%}

%% CHECK INPUT ============================================================

% - Base volume
if ischar(base) % for filenames
  base = load_nifti(base);
end
assert(isstruct(base),'BASE should be a structure!');
if ~isfield(base,'vol') && isfield(base,'.img')  % read from load_untouch_nii
  base.vol = base.img;
  base = rmfield(base,'img');
end
if ~isfield(base,'vox2ras')
  if isfield(base,'hdr')  % read from load_untouch_nii (
    base.vox2ras = [base.hdr.hist.srow_x; ...
      base.hdr.hist.srow_y; base.hdr.hist.srow_z; 0 0 0 1];
  elseif isfield(base,'info') % read from NIFTIINFO (Image Processing Toolbox)
    base.vox2ras = info.Transform.T';
  end
end

% - Data volume
if ischar(data) % for filenames
  data = load_nifti(data);
end
if ~isstruct(data) % numeric vector of matrix
  if isvector(data)
    assert(numel(data) == numel(base.vol), ...
      '1-D data does not seem to be in the same space as base.vol')
    data = reshape(data, size(base.vol));
  end
  data = struct('vol',data, 'vox2ras', base.vox2ras);
end

% OPTION 1
%{
contours for clusters, ROIs, ...
and what other annotations?
%}

%% C O N F I G ============================================================
if ~exist('cfg','var'), cfg=[]; end
if ~isfield(cfg,'method')
  cfg.method = 'linear';
end
if ~isfield(cfg,'basemethod')
  cfg.basemethod = 'linear';
end
if ~isfield(cfg,'contoursmoothing')
  cfg.contoursmoothing = 3;
end

%% COORDINATES
if ~isfield(cfg,'xyz')
  cfg.xyz = 'axi8';
end

% automatic subcortical contours
if strcmpi(cfg.xyz(1:end-1),'subcort')
  nslices = str2double(cfg.xyz(end));
  cfg.xyz = [nan(nslices,1), linspace(-30, 15, nslices)', nan(nslices,1)];

  Labels = [18 54 17 53 26 58 12 51 11 50 13 52 10 49]; % aymgdala, hipppocampus, GP, vStr, ...
  Aseg = helper_conformmri(helper_readmri(fullfile(getenv('SUBJECTS_DIR'),'cvs_avg35_inMNI152','mri','aseg.mgz')));
  cfg.contourcolormap = ones(numel(Labels),3);
  cfg.contourwidth = 1.5;
  cfg.contour = {};
  for iLabel = 1:numel(Labels)
    cfg.contour{iLabel} = Aseg;
    cfg.contour{iLabel}.vol = (cfg.contour{iLabel}.vol == Labels(iLabel));
  end
  CoordsVox = find3(ismember(Aseg.vol, Labels));
  CoordsMm = Aseg.vox2ras*[CoordsVox'-1; ones(1, size(CoordsVox,1))];
  MARGIN_MM = 5;
  cfg.contourbbox = [min(CoordsMm(1:3,:), [], 2)-MARGIN_MM, max(CoordsMm(1:3,:), [], 2)+MARGIN_MM];
end


if ischar(cfg.xyz)
  switch(cfg.xyz(1:3))
    case 'sag'
      xyzdim = 1;
    case 'cor'
      xyzdim = 2;
    case 'axi'
      xyzdim = 3;
  end
  nslices = str2double(cfg.xyz(4:end));
  % equidistance over the DATA volume:
  ijk = [];
  ind_nnz = ~isnan(data.vol) & data.vol~=0;
  [ijk(:,1),ijk(:,2),ijk(:,3)] = ind2sub(size(data.vol), find(ind_nnz));
  ijk = [min(ijk); max(ijk)];
  if ~isempty(ijk)
    bbox = data.vox2ras*[ijk-1 ones(2,1)]';
    xyz = linspace(bbox(xyzdim,1), bbox(xyzdim,2), nslices+2);
    cfg.xyz = nan(nslices,3);
    cfg.xyz(:,xyzdim) = sort(xyz(2:end-1));
  else
    warning('All voxels are masked.')
    cfg.xyz = [nan nan 1];
  end
end
nslices = size(cfg.xyz,1);

%% I'LL HAVE TO THINK MORE ABOUT THE LAYOUT...
% NOW ALWAYS [1 X (#SLICES+1)]:
if cfg.IsHistogram
  cfg.layout = [1 nslices+1];
else
  cfg.layout = [1 nslices];
end
if ~isfield(cfg,'axesinsetmargin')
  cfg.axesinsetmargin = [.1, .03, .07, .13];
end
allaxes = axeslayout(cfg.layout, cfg.axesinsetmargin,'tight');
cfg.sliceaxes = struct('w',allaxes.w, 'h',allaxes.h, 'x',allaxes.x(1:nslices), 'y',allaxes.y(1:nslices));
allaxes = axeslayout(cfg.layout, cfg.axesinsetmargin+[0.2 0 0 0],'tight');
cfg.histaxes = struct('w',allaxes.w, 'h',allaxes.h, 'x',allaxes.x(end), 'y',allaxes.y(end));

%% Figure
if ~isfield(cfg,'figureposition')
  figpos = get(0,'defaultFigurePosition');
  cfg.figureposition = [figpos(1:2)  200*cfg.layout(2) 200];
end

%% Color range
if ~isfield(cfg,'mask')
  cfg.mask = true(size(data.vol));
end
numvals = data.vol(:);
numvals(~cfg.mask(:)) = [];
numvals(isnan(numvals)) = [];
numvals(isinf(numvals)) = [];
if isempty(numvals)
  warning('All vertices are masked')
  numvals = 0;
end
if ~isfield(cfg,'caxis')
  cfg.caxis = winsorcaxis(numvals);
end
if ~isfield(cfg,'thres')
  cfg.thres = [0 0];
end
if numel(cfg.thres) == 1
  cfg.thres = [-abs(cfg.thres) abs(cfg.thres)];
end

%% -- Color map
if ~isfield(cfg,'subthres')
  cfg.subthres = 0;
end
if ~isfield(cfg,'colormap')
  if strcmp(cfg.method,'mip')
    cfg.colormap = flipud(gray);
  else
    cfg.colormap = getcolormap(numvals, cfg);
  end
else % if a colormap is given
  if ~isequal(cfg.thres,[0 0]) % THRESHOLDED?
    cfg.colormap = threscolormap(cfg);
  end
end
if ~isfield(cfg,'basecolormap')
  cfg.basecolormap = gray;
end

%% M A I N ================================================================
%% -- Initialize figure
if ~isfield(cfg,'figurehandle')
  cfg.figurehandle = figure;
else
  if ~isfield(cfg,'figurehold')
    cfg.figurehold = false;
  end
  if cfg.figurehold
    hold on
  else
    clf
  end
end
set(gcf, 'position', cfg.figureposition);
if isfield(cfg,'fname_png') % if fname_png is given, make it invisible
  set(gcf,'visible','off')
end

%% -- DRAW slices
H = struct();

for iaxes = 1:nslices
  % Set axes
  H(iaxes).baseaxes = axespos(cfg.sliceaxes, iaxes);
  H(iaxes).overaxes = axespos(cfg.sliceaxes, iaxes);

  % Get slices
  if size(data.vol,4) == 3
    Vover = [];
    for j = 1:3
      [Vover_, Uover, Wover] = vol2slice(data.vol(:,:,:,j), vox2ras_0to1(data.vox2ras), cfg.xyz(iaxes,:), cfg.method);
      Vover = cat(3, Vover, Vover_);
    end
  else
    [Vover, Uover, Wover] = vol2slice(data.vol(:,:,:,1), vox2ras_0to1(data.vox2ras), cfg.xyz(iaxes,:), cfg.method);
  end
  [Vbase, Ubase, Wbase] = vol2slice(base.vol(:,:,:,1), vox2ras_0to1(base.vox2ras), cfg.xyz(iaxes,:), cfg.basemethod);

  switch (cfg.method)
    case 'mip'
      % - mip of overlay
      H(iaxes).overslice = helper_over(H(iaxes).overaxes, Vover, Uover, Wover, cfg);

      % - base contour
      H(iaxes).baseslice = helper_contour(H(iaxes).baseaxes, Vbase, Ubase, Wbase, cfg);

      % - equalize axes
      helper_equalizeaxes(H(iaxes).baseaxes, H(iaxes).overaxes)

      % - mip-specific setting
      set(H(iaxes).baseaxes,'color','none')
      grid(H(iaxes).baseaxes,'on')

    otherwise
      % - base image
      H(iaxes).baseslice = imagesc(H(iaxes).baseaxes, Ubase.axis, Wbase.axis, Vbase);

      % - overlay
      H(iaxes).overslice = helper_over(H(iaxes).overaxes, Vover, Uover, Wover, cfg);

      % - equalize axes
      helper_equalizeaxes(H(iaxes).baseaxes, H(iaxes).overaxes)

      % - overlay-specific setting
      set(H(iaxes).baseaxes,'xtick',[],'ytick',[])
      grid(H(iaxes).baseaxes,'off')
  end

  % -- WORLD cooridinate (XYZ=RAS)
  xyzdim = find(~isnan(cfg.xyz(iaxes,:)));
  xyzlabel = 'XYZ';
  xlabel(H(iaxes).baseaxes, sprintf('%s = %.0f', xyzlabel(xyzdim), cfg.xyz(iaxes,xyzdim)));


  % - add contours
  if isfield(cfg,'contour')
    if ~isfield(cfg,'contourwidth'), cfg.contourwidth = 1; end
    H(iaxes).contour = [];
    ncons = numel(cfg.contour);
    if ~isfield(cfg,'contourcolormap')
      cfg.contourcolormap = brewermap(ncons, 'Set1');
    end
    for icon = 1:ncons
      [Vi, U, W] = vol2slice(...
        cfg.contour{icon}.vol, vox2ras_0to1(cfg.contour{icon}.vox2ras), cfg.xyz(iaxes,:), 'nearest');
      hold on; warning off
      [~,H(iaxes).contour] = contour(H(iaxes).overaxes, U.axis, W.axis, ...
        convn(Vi,ones(cfg.contoursmoothing),'same'), 1); % slight smoothing
      H(iaxes).contour.Color = cfg.contourcolormap(icon,:);
      H(iaxes).contour.LineWidth = cfg.contourwidth;
      hold off; warning on
    end

    if isfield(cfg, 'contourbbox')
      % - set zoom
      leftdims = setdiff(1:3, xyzdim);
      uIdx = leftdims(1);
      xlim(H(iaxes).overaxes, cfg.contourbbox(uIdx,:))
      xlim(H(iaxes).baseaxes, cfg.contourbbox(uIdx,:))
      wIdx = leftdims(2);
      ylim(H(iaxes).overaxes, cfg.contourbbox(wIdx,:))
      ylim(H(iaxes).baseaxes, cfg.contourbbox(wIdx,:))

    end
  end

  % TODO: - add annotations



  % - common setting
  set(H(iaxes).baseaxes, 'xticklabel',[],'yticklabel',[], 'DataAspectRatio',[1 1 1], 'Ydir','nor')
  set(H(iaxes).overaxes, 'visible','off', 'DataAspectRatio',[1 1 1], 'Ydir','nor')
  colormap(H(iaxes).baseaxes, cfg.basecolormap)
  colormap(H(iaxes).overaxes, cfg.colormap)
  clim(H(iaxes).overaxes, cfg.caxis)

end

%% -- DRAW histogram
if cfg.IsHistogram
  H(iaxes).histaxes = axespos(cfg.histaxes, 1);
  numvals(~numvals) = []; % discard zeros
  H(iaxes).hist = colorhist(numvals, cfg);
  axis square
  if isfield(cfg,'histtitle')
    title(H(iaxes).histaxes, cfg.histtitle, 'interp','none')
  end
end

%% OUT
if isfield(cfg,'fname_png')
  if ~isfield(cfg,'dpi')
    cfg.dpi = 300;
  end
  if ~isempty(getCurrentWorker)
    rendopt = '-painters';
  else
    rendopt = '-opengl';
  end
  export_fig(cfg.fname_png,['-r',num2str(cfg.dpi)],rendopt)
  close(cfg.figurehandle)
end
if ~nargout, clear H cfg; end

end


%% UNIT-TEST

function TEST()
mri = load_nifti('/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz');
thal = load_nifti('/usr/local/fsl/data/atlases/Thalamus/Thalamus-maxprob-thr50-2mm.nii.gz');
prob = load_nifti('/usr/local/fsl/data/atlases/JHU/JHU-ICBM-tracts-prob-2mm.nii.gz');
prob.vol = prob.vol(:,:,:,1) - prob.vol(:,:,:,2);
H = myfs_viewvol(mri, prob, struct('xyz',[nan nan -20; nan nan -10; nan nan 0; nan nan 10; nan nan 20; nan nan 30], ...
  'method','linear','caxis',[-50 50],'contour',{{thal}}, 'contoursmoothing',1) )
%%
surfs = myfs_readsurfs('fsaverage');
mri = load_nifti(fullfile(getenv('FSL_DIR'),'/data/standard/MNI152_T1_1mm.nii.gz'));
prob = load_nifti(fullfile(getenv('FSL_DIR'),'data/atlases/JHU/JHU-ICBM-tracts-prob-2mm.nii.gz'));
prob.vol = prob.vol(:,:,:,1) - prob.vol(:,:,:,2);
spaces = surfs;
spaces.mri = mri;
data = {surfs.thickness{1}-2, surfs.thickness{2}-2, prob};
myfs_viewsurfvol(spaces, data, struct('xyz','axi6','colormap',flipud(brewermap(256,'spectral'))))

end
