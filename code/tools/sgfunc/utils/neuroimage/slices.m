function [H, cfg, base, data] = slices(base, data, cfg)
%
% [H, cfg] = slices(base, data, cfg)
%
% [INPUT]
% BASE can be (1) filename, (2) structure, (3) 3-D numarray
% DATA can be (1) filename, (2) structure, (3) 3-D numarray, (4) vector
%             (5) 'fslho-thr0' | 'fslho-thr25' | 'fslho-thr50'
% CFG is a structure:
% (.contour) can be (1) filename, (2) structure, or (3) 3-D numarray
%                   or a cell array of such
%
% (cc) 2021, dr.seunggoo.kim@gmail.com

%% C O N F I G ============================================================
if ~exist('cfg','var'), cfg=[]; end
if ~isfield(cfg,'method')
  cfg.method = 'linear';
end
if ~isfield(cfg,'basemethod')
  cfg.basemethod = 'linear';
end
if ~isfield(cfg,'contoursmoothing')
  cfg.contoursmoothing = 1;
end

%% CHECK INPUT ============================================================
% - Base volume
if ischar(base) % for filenames
  base = helper_read(base);
end
if islogical(base)
  base = double(base);
end
if isnumeric(base)
  warning('Numeric base volume: vox2ras = eye(4) assumed.')
  Q = eye(4); Q(1:3,4) = 1;
  base = struct('vol',base, 'vox2ras',Q);
end
base = helper_conform(base); % make sure all have .vol and .vox2ras

% - Data volume
if ~exist('data','var'), data = []; end
if ischar(data)
  [p1,f1,e1] = myfileparts(data);
  if contains(lower(data),'fslho') && isempty(e1)
    [data, cmap] = helper_readfslatlas(data);
    if ~isfield(cfg, 'colormap')
      cfg.colormap = cmap;
    end
    cfg.method = 'nearest';
  else
    data = helper_read(data);
  end
end
if islogical(data)
  data = double(data);
end
if isnumeric(data) % numeric vector of matrix
  if ~isempty(data) % Allowing data to be null (possibly only base+contour)
    if isvector(data)
      assert(numel(data)==numel(base.vol), ...
        '1-D data does not seem to be in the same space as base.vol')
      data = reshape(data, size(base.vol));
    end
    data = struct('vol',data, 'vox2ras', base.vox2ras);
  end
end
if ~isempty(data)
  data = helper_conform(data);
end


%% - Contour volume
if isfield(cfg,'contour')
  if ~iscell(cfg.contour)
    cfg.contour = {cfg.contour};
  end
  ncons = numel(cfg.contour);
  for icon = 1:ncons
    if ischar(cfg.contour{icon}) % for filenames
      cfg.contour{icon} = helper_read(cfg.contour{icon});
    end
    if islogical(cfg.contour{icon})
      cfg.contour{icon} = double(cfg.contour{icon});
    end
    if isnumeric(cfg.contour{icon})
      cfg.contour{icon} = struct('vol',cfg.contour{icon}, ...
        'vox2ras',base.vox2ras);
    end
    % make sure all have .vol and .vox2ras:
    cfg.contour{icon} = helper_conform(cfg.contour{icon});
  end
end

%% XYZ coordinates
if ~isfield(cfg,'xyz')
  cfg.xyz = ['sag3';'cor3';'axi3'];
end

% bounding box for BASE:
ijk = [];
ind_nnz = ~isnan(base.vol) & base.vol~=0;
[ijk(:,1),ijk(:,2),ijk(:,3)] = ind2sub(size(base.vol), find(ind_nnz));
ijk = [min(ijk); max(ijk)];
if ~isempty(ijk)
  bbox_base = base.vox2ras*[ijk-1 ones(2,1)]';
else
  bbox_base = [1 1 1 1; size(base.vol) 1]';
end

% bounding box for DATA ???
if ~isempty(data)
  ijk = [];
  ind_nnz = ~isnan(data.vol) & data.vol~=0;
  [ijk(:,1),ijk(:,2),ijk(:,3)] = ind2sub(size(data.vol), find(ind_nnz));
  ijk = [min(ijk); max(ijk)];
  if ~isempty(ijk)
    bbox_data = data.vox2ras*[ijk-1 ones(2,1)]';
  else
    bbox_data = [1 1 1 1; size(data.vol) 1]';
  end
end

% set up coordinates
if ischar(cfg.xyz)
  xyz = [];
  for irow = 1:size(cfg.xyz,1)
    switch(cfg.xyz(irow,1:3))
      case 'sag'
        xyzdim = 1;
      case 'cor'
        xyzdim = 2;
      case 'axi'
        xyzdim = 3;
    end
    nslices = str2double(cfg.xyz(irow,4:end)); % # slices for this row

    % equidistance over a volume (Data if exist; otherwise BASE)
%     if ~isempty(data)
%       bbox = bbox_data;
%     else
%       bbox = bbox_base;
%     end
    bbox = bbox_base; % Why data? when it is useful? when you have sparse data consistently across all conditions, but if not this can be inconvenient.
    coords = linspace(bbox(xyzdim,1), bbox(xyzdim,2), nslices+2);
    xyz_add = nan(nslices,3);
    xyz_add(:,xyzdim) = sort(coords(2:end-1)); % excluding both ends
    xyz = [xyz; xyz_add];
  end
  cfg.xyz = xyz;
end
nslices = size(cfg.xyz,1);

%% ANNOTATIONS
if ~isfield(cfg,'showticks')
  cfg.showticks = false;
end
% background color of the coordinate label, color, size...


%% LAYOUT
if ~isfield(cfg,'layout')
  cfg.layout = [ceil(sqrt(nslices)) ceil(sqrt(nslices))];
end
% if ischar(cfg.layout)
%   cfg.layout = double(strsplit(cfg.layout,'x'));
% end
if nslices > prod(cfg.layout)
  error('nslices > prod(cfg.layout)')
end
if ~cfg.showticks
  cfg.sliceaxes = axeslayout(cfg.layout, [0 0 0 0],[0 0 0 0]);
else
  cfg.sliceaxes = axeslayout(cfg.layout, [0.1 0 0 0.1],[0 0 0 0]);
end

%% Figure
if ~isfield(cfg,'figureposition')
  figpos = get(0,'defaultFigurePosition');
  cfg.figureposition = [figpos(1:2)  150*cfg.layout(2) 150*cfg.layout(1)];
end
if ~isfield(cfg,'figurecolor')
  cfg.figurecolor = 'k';
end

%% Color range
% DATA:
if ~isempty(data)
  if ~isfield(cfg,'mask')
    cfg.mask = true(size(data.vol));
  end
  numvals = data.vol(:);
  numvals(~cfg.mask(:)) = [];
  numvals(numvals==0) = [];
else
  numvals = [];
end
numvals(isnan(numvals)) = [];
numvals(isinf(numvals)) = [];
if isempty(numvals) && ~isempty(data)
  warning('All DATA voxels are masked')
  numvals = 0;
end
if ~isfield(cfg,'caxis')
  [cfg.caxis, ~] = winsorcaxis(numvals);
end
if ischar(cfg.caxis)
  if strcmp(cfg.caxis,'minmax')
    cfg.caxis = [min(numvals) max(numvals)];
  end
end
if ~isfield(cfg,'thres')
  cfg.thres = [0 0];
end
if numel(cfg.thres) == 1
  cfg.thres = [-abs(cfg.thres) abs(cfg.thres)];
end

% BASE:
if ~isfield(cfg,'basecaxis')
  basenumvals = base.vol(:);
  basenumvals(numvals==0) = [];
  [~,cfg.basecaxis] = winsorcaxis(basenumvals);
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
end
if ~isequal(cfg.thres,[0 0])
  cfg.colormap = threscolormap(cfg);
end

%% -- Figure color and basecolormap
if ~isfield(cfg,'figurecolor')
  cfg.figurecolor = 'k';
end
if ~isfield(cfg,'basecolormap')
  cfg.basecolormap = gray;
end

%% -- Annotations
if ~isfield(cfg,'coordfontcolor')
  cfg.coordfontcolor = 'w';
end
if ~isfield(cfg, 'coordfontsize')
  cfg.coordfontsize = 7;
end

%% -- Remove temp variables
clear numvals

%% M A I N ================================================================
%% -- Initialize figure
if ~isfield(cfg,'figurehandle')
  cfg.figurehandle = figure;
else
%   if ~isfield(cfg,'figurehold')
%     cfg.figurehold = false;
%   end
%   if cfg.figurehold
%     hold on
%   else
%     hold off
%   end
end
set(gcf, 'position', cfg.figureposition, 'color', cfg.figurecolor);
if isfield(cfg,'fname_png') % if fname_png is given, make it invisible
  set(gcf,'visible','off')
end

%% -- DRAW slices
H = struct();
for iaxes = 1:nslices
  % Set axes
  H(iaxes).baseaxes = axespos(cfg.sliceaxes, iaxes);
  grid(H(iaxes).baseaxes,'on')

  if ~isempty(data)
    H(iaxes).overaxes = axespos(cfg.sliceaxes, iaxes);
  end

  % Get slices
  [Vbase,Ubase,Wbase] = vol2slice(...
    base.vol, vox2ras_0to1(base.vox2ras), cfg.xyz(iaxes,:), ...
    cfg.basemethod);

  if ~isempty(data)
    [Vover,Uover,Wover] = vol2slice(...
      data.vol, vox2ras_0to1(data.vox2ras), cfg.xyz(iaxes,:), cfg.method);
  end

  % Draw base/over slices
  switch (cfg.method)
    case 'mip'
      % - mip of overlay
      H(iaxes).overslice = helper_over(...
        H(iaxes).overaxes, Vover, Uover, Wover, cfg);

      % - base contour
      H(iaxes).baseslice = helper_contour(...
        H(iaxes).baseaxes, Vbase, Ubase, Wbase, cfg);

      % - mip-specific setting
      set(H(iaxes).baseaxes, 'color','none')
      axis(H(iaxes).baseaxes, 'image')

      % - equalize axes
      equalizeaxes(H(iaxes).baseaxes, H(iaxes).overaxes)

    otherwise
      % - base image
      H(iaxes).baseslice = helper_img(...
        H(iaxes).baseaxes, Vbase, Ubase, Wbase);

      % - overlay
      if ~isempty(data)
        H(iaxes).overslice = helper_over(...
          H(iaxes).overaxes, Vover, Uover, Wover, cfg);
        equalizeaxes(H(iaxes).baseaxes, H(iaxes).overaxes)
      end

      % - overlay-specific setting
      if ~cfg.showticks
        set(H(iaxes).baseaxes,'xtick',[],'ytick',[])
        grid(H(iaxes).baseaxes,'off')
      end

  end

  % - add contours
  if isfield(cfg,'contour')
    H(iaxes).contour = {};
    ncons = numel(cfg.contour);
    if ~isfield(cfg,'contourwidth'), cfg.contourwidth = 1; end
    if isfield(cfg,'contourcolormap')
      cmap = cfg.contourcolormap;
    else
      cmap = brewermap(ncons, 'Set1');
    end
    for icon = 1:ncons
      [Vi,U,W] = vol2slice(...
        cfg.contour{icon}.vol, vox2ras_0to1(cfg.contour{icon}.vox2ras), ...
        cfg.xyz(iaxes,:), 'nearest');
      if ~isempty(data)
        axes_ref = H(iaxes).overaxes;
      else
        axes_ref = H(iaxes).baseaxes;
      end
      H(iaxes).contour{icon} = helper_contour(axes_ref, Vi, U, W, cfg);
      H(iaxes).contour{icon}.Color = cmap(icon,:);
      H(iaxes).contour{icon}.LineWidth = cfg.contourwidth;
    end
  end

  % - add annotations

  % -- WORLD cooridinate (XYZ=RAS)
  xyzdim = find(~isnan(cfg.xyz(iaxes,:)));
  xyzlabel = 'XYZ';
  u1 = median(Ubase.axis);
  w1 = prctile(Wbase.axis,95);
  text(H(iaxes).baseaxes, u1, w1, ...
    sprintf('%s = %.0f', xyzlabel(xyzdim), cfg.xyz(iaxes,xyzdim)), ...
    'fontsize',cfg.coordfontsize, 'color',cfg.coordfontcolor, ...
    'HorizontalAlignment','center');

  % -- FILE NAME? TITLE BAR?


  % - common setting
  set(H(iaxes).baseaxes, ...
    'DataAspectRatio',[1 1 1], 'Ydir','nor', 'Visible','off')
  colormap(H(iaxes).baseaxes, cfg.basecolormap)
  if cfg.showticks
    xlabel(H(iaxes).baseaxes, Ubase.axisname);
    ylabel(H(iaxes).baseaxes, Wbase.axisname);
    H(iaxes).baseaxes.FontSize = 7;
  end

  if ~isempty(data)
    set(H(iaxes).overaxes, 'visible','off', ...
      'DataAspectRatio',[1 1 1], 'Ydir','nor')
    colormap(H(iaxes).overaxes, cfg.colormap)
    try
      caxis(H(iaxes).overaxes, cfg.caxis)
    catch ME
      warning('caxis not sane')
    end
    try
      caxis(H(iaxes).baseaxes, cfg.basecaxis)
    catch ME
      warning('caxis not sane')
    end
  else
    try
      caxis(H(iaxes).baseaxes, cfg.basecaxis)
    catch ME
      warning('caxis not sane')
    end
  end
end


%% Colorbar
iaxes = min(prod(cfg.layout), ceil(nslices/2)+1);
if ~isempty(data)
  apos = get(H(iaxes).overaxes,'position');
  H(iaxes).colorbar = colorbar(H(iaxes).overaxes);
  set(H(iaxes).overaxes,'position',apos)
else
  apos = get(H(iaxes).baseaxes,'position');
  H(iaxes).colorbar = colorbar(H(iaxes).baseaxes);
  set(H(iaxes).baseaxes,'position',apos)
end
H(iaxes).colorbar.Color = cfg.coordfontcolor;
H(iaxes).colorbar.Location = 'southOutside';
if isfield(cfg,'colorbarposition')
  H(iaxes).colorbar.Position = cfg.colorbarposition;
else
  H(iaxes).colorbar.Position = [.4 .41 .2 .01];
end
  
H(iaxes).colorbar.FontSize = cfg.coordfontsize;
if isfield(cfg,'colorbarxlabel')
  H(iaxes).colorbar.Label.String = cfg.colorbarxlabel;
end
if isfield(cfg,'colorbartitle')
  title(H(iaxes).colorbar, cfg.colorbartitle, 'color',cfg.coordfontcolor)
end
if isfield(cfg,'colorbarvisible')
  H(iaxes).colorbar.Visible = cfg.colorbarvisible;
end

%% Title
if isfield(cfg,'title')
%   title(H(iaxes).colorbar, cfg.title, 'color', cfg.coordfontcolor)
  title(cfg.title, 'color', cfg.coordfontcolor)
end

%% OUT
if isfield(cfg,'fname_png')
  if ~isfield(cfg,'dpi')
    cfg.dpi = 300;
  end
  gl = opengl('data');
  if strcmp(gl.Renderer, 'None')
    rendopt = '-painters';
  else
    rendopt = '-opengl';
  end
  export_fig(cfg.fname_png,['-r',num2str(cfg.dpi)],rendopt)
  close(cfg.figurehandle)
end
if ~nargout, clear H cfg; end

end

%-------------------------------------------------------------------------%
% SUBFUNCTIONS                                                            %
%-------------------------------------------------------------------------%

function mri = helper_read(fname)
[~,~,ext] = fileparts_gz(fname);
switch (ext)
  case {'.nii','.nii.gz','.img'}
    % Try "my" version of MATLAB Image Processing Toolbox (since 2016)
    if exist('niftiinfogz','file') && exist('niftireadgz','file')
      [V,info] = niftireadgz(fname);
      mri = struct('vol',V, 'info',info);
      % NOTE: NIFTIREAD doesn't upscale precision (hmm...)

    elseif exist('load_nifti','file') % FreeSurfer
      mri = load_nifti(fname);
    elseif exist('load_untouch_nii','file') % NIFTI toolbox
      mri = load_untouch_nii(fname);
    else
      error('CANNOT FIND any function to read NIFTI/ANALYZE files!')
    end
  case {'.mgh','.mgz'}
    [~, M, P] = load_mgh(fname, [], [], 1);
    mri = struct('vox2ras',M, 'tr',P(1), 'filpangle',P(2), 'te',P(3), ...
      'ti',P(3), 'fov',P(4));
    [mri.vol] = load_mgh(fname);
  otherwise
    error('EXTENSION UNRECOGNIZED: %s',ext)
end
mri.vol = double(mri.vol); % now I just upscale it as I won't save it
end

function mri = helper_conform(mri)
if ~isfield(mri,'vol')
  if isfield(mri,'img')  % from LOAD_UNTOUCH_NII
    mri.vol = mri.img;
    mri = rmfield(mri,'img');
  else
    error('UNKNOWN volume fieldname??')
  end
end
if ~isfield(mri,'vox2ras')
  if isfield(mri,'hdr')  % from LOAD_UNTOUCH_NII
    mri.vox2ras = [mri.hdr.hist.srow_x; ...
      mri.hdr.hist.srow_y; mri.hdr.hist.srow_z; 0 0 0 1];
  elseif isfield(mri,'info')  % from NIFTIINFO
    mri.vox2ras = mri.info.Transform.T';  % 0-based
  else
    error('UNKNOWN header fieldname??')
  end
end
if size(mri.vol,4) > 1
  fprintf('[%s] 4-D image is given. Showing the only first volume.\n', ...
    mfilename);
  mri.vol = mri.vol(:,:,:,1);
end
mri.vol = double(mri.vol);

end

function h = helper_over(Ha, Vi, u, w, cfg)
if isequal(cfg.thres, [0 0]) % no thresholding
  Vi(Vi==0) = nan;
  h = pcolor(Ha, u.axis, w.axis, Vi);
  h.LineStyle = 'none';
  if isfield(cfg,'facealpha')
    h.FaceAlpha = cfg.facealpha;
  end
else
  % suprathreshold image:
  Vsupra = Vi;
  Vsupra(Vsupra==0) = nan;
  Vsupra(cfg.thres(1)<Vi & Vi<cfg.thres(2)) = nan;
  h = pcolor(Ha, u.axis, w.axis, Vsupra);
  h.LineStyle = 'none';

  % subthreshold image:
  if cfg.subthres
    hold on
    Vsub = Vi;
    Vsub(Vsub==0) = nan;
    Vsub(Vi<cfg.thres(1) | cfg.thres(2)<Vi) = nan;
    h2 = pcolor(Ha, u.axis, w.axis, Vsub);
    h2.LineStyle = 'none';
    h2.FaceAlpha = 0.5;
    hold off
  end
end
set(Ha,'color','none')
end

function Hi = helper_img(Ha, Vi, u, w)
Hi = imagesc(Ha, u.axis, w.axis, Vi);
end

function Hc = helper_contour(Ha, Vi, u, w, cfg)
hold on
if cfg.contoursmoothing % slight smoothing?
[~,Hc] = contour(Ha, u.axis, w.axis, ...
  convn(Vi,ones(cfg.contoursmoothing),'same'), 1); 
else % no smoothing
  [~,Hc] = contour(Ha, u.axis, w.axis, Vi, 1);
end
hold off
end

function helper_annot()

end

function equalizeaxes(baseaxes, overaxes)
xlim = [min(baseaxes.XLim(1),overaxes.XLim(1)), ...
  max(baseaxes.XLim(2),overaxes.XLim(2))];
ylim = [min(baseaxes.YLim(1),overaxes.YLim(1)), ...
  max(baseaxes.YLim(2),overaxes.YLim(2))];
baseaxes.XLim = xlim;
baseaxes.YLim = ylim;
overaxes.XLim = xlim;
overaxes.YLim = ylim;
end

function [data, cmap] = helper_readfslatlas(atlasdescp)
% Blame windows users:
if ~isunix && ~ismac, error('Why Windows?'), end 

%--- FSL Harvard-Oxford Cort+Subcort --------------------------------------

% Find filenames
fslpath = getenv('FSLDIR');
str = strsplit(atlasdescp,'-');
if numel(str)>1
  suffix = str{2};
else
  suffix = 'thr25';
end
fn_atl{1} = [fslpath,'/data/atlases/HarvardOxford/',...
  'HarvardOxford-cort-maxprob-',suffix,'-1mm.nii.gz'];
fn_xml{1} = [fslpath,'/data/atlases/HarvardOxford-Cortical.xml'];
fn_atl{2} = [fslpath,'/data/atlases/HarvardOxford/',...
  'HarvardOxford-sub-maxprob-',suffix,'-1mm.nii.gz'];
fn_xml{2} = [fslpath,'/data/atlases/HarvardOxford-Subcortical.xml'];
fn = dir([fslpath,'/fslpython/envs/fslpython/lib/python*/',...
  'site-packages/fsleyes/assets/luts/harvard-oxford-cortical.lut']);
assert(numel(fn)==1)
fn_lut{1} = fullfile(fn.folder, fn.name);
fn = dir([fslpath,'/fslpython/envs/fslpython/lib/python*/',...
  'site-packages/fsleyes/assets/luts/harvard-oxford-subcortical.lut']);
assert(numel(fn)==1)
fn_lut{2} = fullfile(fn.folder, fn.name);
for i = 1:2
  assert(isfile(fn_atl{i}), 'file "%s" NOT FOUND!', fn_atl{i})
  assert(isfile(fn_xml{i}), 'file "%s" NOT FOUND!', fn_xml{i})
  assert(isfile(fn_lut{i}), 'file "%s" NOT FOUND!', fn_lut{i})
end

% Read files
nctx = 48;
nstx = 21;
ctx = helper_read(fn_atl{1});
assert(max(ctx.vol(:))==nctx)
stx = helper_read(fn_atl{2});
assert(max(stx.vol(:))==nstx)
% remove large masks (cortical-WM, cortcial-GM, CSF)
lbl2remove = [1 2 3 12 13 14];
stx.vol(ismember(stx.vol, lbl2remove)) = 0;

% Combine CORTICAL + SUBCORTICAL
data = ctx;
data.vol(stx.vol>0) = nctx+stx.vol(stx.vol>0);
data.vol(data.vol==0) = nan;

% Read lookup table:
tbl_ctx = myfsl_readlut(fn_lut{1});
assert(size(tbl_ctx,1) == nctx)
tbl_stx = myfsl_readlut(fn_lut{2});
assert(size(tbl_stx,1) == nstx)

cmap_ctx = [tbl_ctx.r, tbl_ctx.g, tbl_ctx.b];
cmap_stx = [tbl_stx.r, tbl_stx.g, tbl_stx.b];
cmap = [cmap_ctx; cmap_stx];



end
