function [H, cfg] = myfs_viewsurfvol (spaces, data, cfg)
%MYFS_VIEWSURFVOL visualizes scalar data on cortical surfaces and volume.
%

%% C O N F I G U R E A T I O N ============================================
if ~exist('cfg','var')
  cfg = [];
end

% -- Basesurf
if ~isfield(cfg,'basesurf'), cfg.basesurf = 'inflated'; end

%% Read volume data for given filenames:
if ischar(spaces.mri)
   spaces.mri = helper_read(spaces.mri);
end
if ischar(data{3})
   data{3} = helper_read(data{3});
end

%% Conform volume data:
spaces.mri = helper_conform(spaces.mri);
data{3} = helper_conform(data{3});

%% CHECK input
nverts = [size(spaces.(cfg.basesurf){1}.vertices,1), ...
    size(spaces.(cfg.basesurf){2}.vertices,1) numel(spaces.mri.vol)];
if isempty(data)
  warning('data is null')
  data = {nan(nverts(1),1), nan(nverts(2),1), nan(size(spaces.mri.vol))};
end
if ~iscell(data) && isvector(data)
  assert(numel(data)==sum(nverts),'1-D data is inconsistent with spaces');
  data = {data(1:nverts(1)), data(nverts(1)+1:nverts(1)+nverts(2)), ...
    reshape(data(sum(nverts(1:2))+1:end), size(spaces.mri.vol))};
end

%% CREATE COLORMAP & CAXIS BASED on SURFACE & VOLUME:
numvals = [data{1}(:); data{2}(:); data{3}.vol(:)];
if isfield(cfg,'masks')
  numvals(~[cfg.masks{:}]) =[];
end
if isfield(cfg,'mask')
  numvals(~cfg.mask(:)) =[];
end
if ~isfield(cfg,'caxis')
  cfg.caxis = winsorcaxis(numvals);
end


%% DRAW surfaces
cfg1 = cfg;
cfg1.layout = '1x4';
if isfield(cfg,'fname_png')
  cfg1 = rmfield(cfg1, 'fname_png');
end
cfg1.figurehandle = figure('visible','off');
[h1, cfg1] = myfs_viewsurf(spaces, data(1:2), cfg1);


%% DRAW slices
% give me some space:
set(cfg1.figurehandle,'position', get(gcf,'position').*[1 1 1 1.6])
for i = 1:4
  h1(i).axes.Position = h1(i).axes.Position+[0 +0.15 0 0];
end
if numel(h1) > 4
  for i = 5:numel(h1) % histograms & colorbar?
    h1(i).axes.Position = h1(i).axes.Position+[0 0.35 0 -0.05];
  end
end

cfg2 = cfg1;
cfg2.basecolormap = gray;
cfg2.overcolormap = cfg1.colormap;
cfg2.figurehandle = gcf;
cfg2.figureposition = get(gcf,'position');
cfg2.figurehold = true;
if ~isfield(cfg2,'xyz')
  cfg2.xyz = 'axi8';
end
if isfield(cfg, 'volaxesinsetmargin')
  cfg2.axesinsetmargin = cfg.volaxesinsetmargin;
else
cfg2.axesinsetmargin = [.05, .05, .60, .05];
end

[h2, cfg2] = myfs_viewvol(spaces.mri, data{3}, cfg2);


%% OUTPUT
if isfield(cfg,'fname_png')
  if ~isfield(cfg,'dpi'), cfg.dpi = 300; end
  export_fig(cfg2.figurehandle, cfg.fname_png, cfg.dpi)
  close(cfg2.figurehandle)
else
  cfg2.figurehandle.Visible = 'on';
end

if nargout > 0
  H = {h1, h2};
end
if nargout > 1
  cfg = cfg2;
end


end



%% In-function functions


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
    mri = struct('vox2ras',M, 'tr',P(1), 'filpangle',P(2), 'te',P(3), 'ti',P(3), 'fov',P(4));
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
    mri.vox2ras = [mri.hdr.hist.srow_x; mri.hdr.hist.srow_y; mri.hdr.hist.srow_z; 0 0 0 1];
  elseif isfield(mri,'info')  % from NIFTIINFO
    mri.vox2ras = mri.info.Transform.T';  % 0-based
  else
    error('UNKNOWN header fieldname??')
  end
end
% if size(mri.vol,4) > 1
%   fprintf('[%s] 4-D image is given. Showing the only first volume.\n', mfilename);
%   mri.vol = mri.vol(:,:,:,1);
% end
mri.vol = double(mri.vol);
end



function TEST()
surfs = myfs_readsurfs('fsaverage4')
myfs_viewsurf(surfs,[])
myfs_viewsurfvol(surfs,randn(2562+2562+182*218*182,1))
myfs_viewsurfvol(surfs,randn(2562+2562+182*218*182,1)+100)
end
