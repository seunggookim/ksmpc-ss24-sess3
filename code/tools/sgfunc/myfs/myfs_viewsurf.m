function [H, cfg] = myfs_viewsurf (surfs, data, cfg)
%MYFS_VIEWSURF visualizes scalar data on cortical surfaces
%
% [USAGE]
% myfs_viewsurf (surfs, data)
% H = myfs_viewsurf (surfs, data)
% [H, cfg] = myfs_viewsurf (surfs, data, cfg)
%
%
% [INPUT]
% surfs {1x2} freesurfer surfaces structure read by MYFS_READSURFS
%
% data  {1x2} cells contain a vector (vertex-mapped scalar data) for each
%             hemisphere. [NaN] will be displayed as transparent if openGL
%             is available.
%     | [Nx1] a vector from the left and right (in this order) hemispheres
%
% cfg   (1x1) structure to configure options (optional):
% - figure
% .figurehandle [1x1] default = creating a new figure
%
% - surfaces layout
%  .basesurf '1xN'  'inflated' (default) | 'smoothwm' | 'white' | 'pial'
%                   ('semiinflated' | 'sphere' | 'suptemp' )
%  .layout   '1xN'  '2x4' (default) | '1x4' | '2x2' | '1x2' | '1x2oblq' |
%                   '1x2stmp'
%  .views    {1xN}  view angle, default is determined by layout:
%                   '2x4' {[-90 0],[90 0],[-90 0],[90 0],
%                          [-90 75],[90 -75],[-90 -75],[90 75]}
%                   '1x4' {[-90 0],[90 0],[-90 0],[90 0]}
%                   '2x2' {[-90 0],[90 0],[90 0],[-90 0]}
%                   '1x2' {[-90 0],[90 0]}
%                   '1x2oblq' {[-119 6],[106 -5]}
%                   '1x2stmp' {[0 90],[0 90]};
%  .hemis    [1xN]  1=left 2=right, default is determined by layout:
%                   '2x4' [1 1 2 2 1 1 2 2 ]
%                   '1x4' [1 1 2 2 ]
%                   '2x2' [1 2 1 2]
%                   '1x2' [1 2]
%  .figureposition  [1x4] default is determined by layout:
%                   '2x4' [5   694   900   365]
%                   '1x4' [5   694   900   220]
%                   '2x2' [5   694   900   793]
%                   '1x2' [5   694   900   415]
%                   '1x2stmp' [5   694   900   615]
%  .figurecolor [1x3] | '1x1' default = 'k' if .thres ~= 0
%                                     = 'w' if .thres == 0
%
% - color schemes (presets)
%  .colorscheme    '1xN'  'darkspectral' | 'darkparula'
% - color coding (manual)
%  .basecolor '1xN'  'dark' (0.15/0.35) | 'bright' (0.7/1)
%  .facealpha [1x1]  default = 1 (overlay)
%  .colormap  [Nx3]  default = [bluish; gray; redish] if .thres ~= 0;
%                    brewermap(256,'spectral') if .thres == 0
%
% -- optional colormaps:
%  .caxis     [1x2] default is determined by signs of values:
%                   for positive values: [1%tile 99%tile]
%                   for sigend values: [-/+max(99%tile of abs)]
%             alternatively, enter 'full' min/max of numerical values
%  .thres     [1x1] | [1x2]  default=0 (no thresholding)
%  .subthres  [1x1]  true | {false} - show subthreshold values with a=0.5
%  .masks     {1x2} cells contain vertex-mapped logical vectors for masking
%                   by default, non-cortical vertices are masked.
%
% - contours
% -- isocurvature contours:
%  .doisocurv        [1x1]  true (requires surf.(basesurf).isocurv)
%                           | false (default)
%  .isocurvcolor     [1x3]  default = [0 0 0]
%  .isocurvlinewidth [1x1]  default = 0.5
%  .isocurvlinestyle '1xN'  default = '-'
% -- cluster contours:
%  .doclusid         [1x1]  true (requires surf.(basesurf).isoclus)
%                           | false (default)
%  .isocluscolor     [1x3]  default = [1 1 1]
%  .isocluslinewidth [1x1]  default = 2
%  .isocluslinestyle '1xN'  default = '-'
% -- outerboundary contours:
%  .outerboundary      [1x1]  true | false (default)
%  .outerboundarycolor [1x3]  default = [0 0 0]
%
% - lighting
%  .camlight    [1x1]  {true} | false
%  .camposition [1x2]  default=[0 10]
%
% - histograms
%  .histfontsize  [1x1]  default=8
%  .histfontsize
%  .histxtick
%  .histxticklabel
%  .histxlim
%
% - colorbar
%  .colorbartitle
%  .colorbarinterp   '1xN'  {'none'} | 'tex' | 'latex'
%  .colorbarfontsize
%  .colorbarxlabel
%  .colorbarxlabelinterp '1xN'  {'none'} | 'tex' | 'latex'
%
% - print (requires EXPORT_FIG)
%  .fname_png
%  .dpi
%
%
% [OUTPUT]
% H structure array (1xN) contains axes handles for:
% - surface axes
%  .axes
%  .basesurf
%  .oversurf
%  .light
%
% - histogram axes
%  .axes
%  .histfontsize
%  .histxtick
%  .histxticklabel
%
% - colorbar axes
%  .axes
%  .colorbar
%  .colorbarxtick
%  .colorbarxticklabel
%  .colorbartitle
%  .colorbarfontsize
%
% (cc) 2019-2020, sgKIM. mailto://solleo@gmail.com
% This is distributed from https://github.com/solleo/surfviz
%
% SEE ALSO: VIEW_SURFDATA, MYFS_READSURFS, EXPORT_FIG, BREWERMAP

%{
Creative Commons Legal Code

CC0 1.0 Universal

    CREATIVE COMMONS CORPORATION IS NOT A LAW FIRM AND DOES NOT PROVIDE
    LEGAL SERVICES. DISTRIBUTION OF THIS DOCUMENT DOES NOT CREATE AN
    ATTORNEY-CLIENT RELATIONSHIP. CREATIVE COMMONS PROVIDES THIS
    INFORMATION ON AN "AS-IS" BASIS. CREATIVE COMMONS MAKES NO WARRANTIES
    REGARDING THE USE OF THIS DOCUMENT OR THE INFORMATION OR WORKS
    PROVIDED HEREUNDER, AND DISCLAIMS LIABILITY FOR DAMAGES RESULTING FROM
    THE USE OF THIS DOCUMENT OR THE INFORMATION OR WORKS PROVIDED
    HEREUNDER.
%}

%{
:HISTORY:
2019-12-17: removed "useparula" option from FSSS_VIEW and VIEW_SURFDATA
2020-02-23: new color schemes: "brightspectral", "brightparula",
"yellowblue"
2020-03-24: isocontour:edge-connectig:0.07 sec/hemi!
2020-06-03: outerboundary
2020-11-12: rename FSSS to MYFS and repackaging
%}
DEBUG = false;

%% P A T H ================================================================
myfs_addpath()

%% C O N F I G ============================================================
if ~exist('cfg','var'), cfg=[]; end


%% -- Basesurf
if ~isfield(cfg,'basesurf')
  if isfield(surfs,'basesurf')
    cfg.basesurf = surfs.basesurf;
  else
    fields = fieldnames(surfs);
    idx = find(ismember(fields, {'semiinflated','inflated','white','pial'}));
    if isempty(idx)
      error('No BASESURF recognized/defined')
    end
    cfg.basesurf = fields{idx(1)};
  end
end


%% -- Input dimension
nverts = [size(surfs.(cfg.basesurf){1}.vertices,1), ...
  size(surfs.(cfg.basesurf){2}.vertices,1)];

% Handle null input:
if isempty(data)
  warning('data is null')
  data = {nan(nverts(1),1), nan(nverts(2),1)};
end

% Handle concatenated vector:
if ~iscell(data)
  if length(data) == sum(nverts)
    warning('data is assumed to be bi-hemi.')
    data = {reshape(data(1:nverts(1)),[],1), ...
      reshape(data((1:nverts(2))+nverts(1)),[],1)};
  end
end

% Forcing an input vector as a column vector:
for ihemi = 1:2
  if isvector(data{ihemi})
    data{ihemi} = reshape(data{ihemi},[],1);
  else
    if size(data{ihemi},2) ~=3
      error('Matrix input should be RGB (3-dimensional)')
    end
  end
end


%% -- SUPTEMP
if contains(cfg.basesurf,'suptemp')
  if ~isfield(cfg,'layout')
    cfg.layout = '1x2stmp';
  end
  if ~isfield(cfg,'masks')
    cfg.masks = {surfs.aparc{1}.suptemp1_bin', surfs.aparc{2}.suptemp1_bin'};
  end
end

%% -- FLAT surfaces
if contains(cfg.basesurf,'_flat')
  if ~isfield(cfg,'layout')
    cfg.layout = '1x2flat';
  end
end

%% -- Color scheme
if isfield(cfg,'colorscheme')
  switch cfg.colorscheme
    case 'darkspectral'
      cfg.colormap = flipud(brewermap(256,'spectral'));
      cfg.figurecolor = 'k';
      cfg.basecolor = 'dark';
    case 'darkparula'
      cfg.colormap = parula(256);
      cfg.figurecolor = 'k';
      cfg.basecolor = 'dark';
    case 'brightspectral'
      cfg.colormap = flipud(brewermap(256,'spectral'));
      cfg.figurecolor = 'w';
      cfg.basecolor = 'bright';
    case 'brightparula'
      cfg.colormap = parula(256);
      cfg.figurecolor = 'w';
      cfg.basecolor = 'bright';
    case 'yellowblue'
      cfg.colormap = [
        [linspace(0,0,64); linspace(1,0.25,64); linspace(1,1,64)]';
        [linspace(1,1,64); linspace(0.25,1,64); linspace(0,0,64)]'];
      cfg.figurecolor = 'k';
      cfg.basecolor = 'dark';
    otherwise
      warning('%s is not valid colorscheme. Ignored', cfg.colorscheme)
  end
end

%% -- Layout
if ~isfield(cfg,'IsColorbar'), cfg.IsColorbar = 1; end
if ~isfield(cfg,'IsHistogram'), cfg.IsHistogram = 1; end

if ~isfield(cfg,'layout'), cfg.layout='2x4'; end
if ~isfield(cfg,'views')
  switch cfg.layout
    case '2x4'
      cfg.views = {...
        [-90 0],[90 0],[-90 0],[90 0], ...
        [-90 75],[90 -75],[-90 -75],[90 75]};
    case '1x4'
      cfg.views = {[-90 0],[90 0],[-90 0],[90 0]};
    case '2x2'
      cfg.views = {[-90 0],[90 0],[90 0],[-90 0]};
    case '3x2'
      cfg.views = {[-90 0],[90 0],[90 0],[-90 0],[180 0],[180 0]};
    case '2x2oblqdn'
      cfg.views = {[-120 10],[120 10],[120 10],[-120 10]};
    case 'audhip'
      cfg.views = {[-90 25],[90 25],[90 -45],[-90 -45]};
    case {'1x2big','1x2','2x1'}
      cfg.views = {[-90 0],[90 0]};
    case '1x2oblq'
      cfg.views = {[-120 10],[120 10]};
    case {'1x2stmp','1x2flat'}
      cfg.views = {[0 90],[0 90]};
  end
end
if ~isfield(cfg,'hemis')
  switch cfg.layout
    case '2x4'
      cfg.hemis = [1 1 2 2 1 1 2 2];
    case '1x4'
      cfg.hemis = [1 1 2 2];
    case {'2x2','2x2oblqdn','audhip'}
      cfg.hemis = [1 2 1 2];
    case '3x2'
      cfg.hemis = [1 2 1 2 1 2];
    otherwise
      if contains(cfg.layout, '1x2','2x1')
        cfg.hemis = [1 2];
      end
  end
end
if ~isfield(cfg,'axes')
  switch cfg.layout
    case '2x4'
      cfg.surfaxes = axeslayout([2 4],[.02 .02 0 0],[.02, .02, .0, .0]);
      cfg.surfaxes.y(1:4) = cfg.surfaxes.y(1:4)+0.02;
      cfg.surfaxes.y(5:8) = cfg.surfaxes.y(5:8)+0.15;
      cfg.histaxes = {[.05 .07 .25 .15],[.7 .07 .25 .15]};
      cfg.colorbaraxes = [.5-.15 .09 .3 .08];
    case '1x4'
      cfg.surfaxes = axeslayout([1 4],[.02 .02 0 0],[.02, .02, .0, .0]);
      cfg.surfaxes.y(1:4) = cfg.surfaxes.y(1:4)+0.15;
      cfg.histaxes = {[.05 .13 .25 .15],[.7 .13 .25 .15]};
      cfg.colorbaraxes = [.5-.15 .15 .3 .08];
    case {'2x2'}
      cfg.surfaxes = axeslayout([2 2],[.02 .02 0 0],[.02, .02, .0, .0]);
      cfg.surfaxes.y(1:2) = cfg.surfaxes.y(1:2)+0.05;
      cfg.surfaxes.y(3:4) = cfg.surfaxes.y(3:4)+0.15;
      cfg.histaxes = {[.07 .04 .25 .15],[.72 .04 .25 .15]};
      cfg.colorbaraxes = [.5-.15 .06 .3 .08];
    case '2x2oblqdn'
      cfg.surfaxes = axeslayout([2 2],[-.08 -.08 .0 .0],[.02, .02, .0, .0]);
      cfg.surfaxes.y(1:2) = cfg.surfaxes.y(1:2)+0.03;
      cfg.surfaxes.y(3:4) = cfg.surfaxes.y(3:4)+0.11;
      cfg.surfaxes.x([1 3]) = cfg.surfaxes.x([1 3])+0.02;
      cfg.surfaxes.x([2 4]) = cfg.surfaxes.x([2 4])-0.02;
      cfg.surfaxes.h = cfg.surfaxes.h*1.03;
      cfg.histaxes = {[.10 .04 .25 .1],[.7 .04 .25 .1]};
      cfg.colorbaraxes = [.5-.15 .06 .3 .08];
    case '3x2'
      cfg.surfaxes = axeslayout([3 2],[.02 .02 0 0],[.02, .02, .0, .0]);
      cfg.histaxes = {[.05 .04 .25 .15],[.7 .04 .25 .15]};
      cfg.colorbaraxes = [.5-.15 .06 .3 .08];
    case {'1x2'}
      cfg.surfaxes = axeslayout([1 2],[.02 .02 0 0],[.02, .02, .0, .2]);
      cfg.surfaxes.y(1:2) = cfg.surfaxes.y(1:2)+0.06;
      cfg.histaxes = {[.08 .08 .23 .15],[.73 .08 .23 .15]};
      cfg.colorbaraxes = [.5-.15 .20 .3 .05];
    case {'1x2big','1x2oblq','1x2flat'}
      cfg.surfaxes = axeslayout([1 2],[.02 .02 0 0],[.02, .02, .0, .0]);
      cfg.surfaxes.y = cfg.surfaxes.y + 0.11;
      cfg.histaxes = {[.05 .07 .25 .15],[.7 .07 .25 .15]};
      cfg.colorbaraxes = [.5-.15 .09 .25 .05];
    case '1x2stmp'
      cfg.surfaxes = axeslayout([1 2],[.02 .02 0 0],[.02, .02, .0, .0]);
      cfg.surfaxes.y = cfg.surfaxes.y + 0.11;
      cfg.histaxes = {[.1 .05 .22 .10],[0.73 .05 .22 .10]};
      cfg.colorbaraxes = [.5-.25/2 .04 .25 .08];
    case 'audhip'
      cfg.surfaxes = axeslayout([2 2],[.02 .02 0 0],[.02, .02, .0, .0]);
      cfg.surfaxes.y(1:2) = cfg.surfaxes.y(1:2)+0.04;
      cfg.surfaxes.y(3:4) = cfg.surfaxes.y(3:4)+0.16;
      cfg.histaxes = {[.05 .06 .25 .15],[.7 .06 .25 .15]};
      cfg.colorbaraxes = [.35 .08 .3 .08];
    case '2x1'
      cfg.surfaxes = axeslayout([2 1],[.02 .02 0 0],[.02, .02, .1, .25]);
      cfg.surfaxes.y(1:2) = cfg.surfaxes.y(1:2)+0.06;
      cfg.histaxes = {[.08 .08 .23 .15],[.73 .08 .23 .15]};
      cfg.colorbaraxes = [.5-.15 .20 .3 .05];
  end
end
if ~isfield(cfg,'figureposition')
  figpos = get(0,'defaultFigurePosition');
  switch cfg.layout
    case '2x4'
      cfg.figureposition = [figpos(1:2)   900   365];
    case '1x4'
      cfg.figureposition = [figpos(1:2)   900   220];
    case '1x2'
      cfg.figureposition = [figpos(1:2)   450   220];
    case {'2x2','2x2oblqdn'}
      cfg.figureposition = [figpos(1:2)   450   365];
    case '3x2'
      cfg.figureposition = [figpos(1:2)   450   600];
    case {'audhip'}
      cfg.figureposition = [figpos(1:2)   450   400];
    case {'1x2big','1x2oblq','1x2flat'}
      cfg.figureposition = [figpos(1:2)   900   415];
    case '1x2stmp'
      cfg.figureposition = [figpos(1:2)   250   360];
    case '2x1'
      cfg.figureposition = [figpos(1:2)   250   500];
  end
end


%% -- Color range based on both hemisphere
if ~isfield(cfg,'masks')
  for hemi = 1:2
    cfg.masks{hemi} = true(size(data{hemi}));
  end
end
for hemi = 1:2
  cfg.masks{hemi} = cfg.masks{hemi} & surfs.aparc{hemi}.cortex;
end
if islogical(data{1,1}) || islogical(data{1,1})
  data = {uint8(data{1}), uint8(data{2})};
end

% get numeric values:
numvals = cat(1,data{:});
numvals(~cat(1,cfg.masks{:})) = [];
numvals(isnan(numvals)) = [];
numvals(isinf(numvals)) = [];
if isempty(numvals)
  warning('No numeric values are given!')
  numvals = 0;
end
if ~isfield(cfg,'isinterger')
  cfg.isinteger = isequal(numvals, round(numvals));
end

% find (mirroed) Winsorized color limits for continous values (not for
% integers):
if ~isfield(cfg,'caxis')
  cfg.caxis = winsorcaxis(numvals);
end

% Define thresholds (or not);
if ~isfield(cfg,'thres')
  cfg.thres = [0 0];
end
if numel(cfg.thres) == 1
  cfg.thres = [-abs(cfg.thres) abs(cfg.thres)];
end

%% -- Colormap
if ~isfield(cfg,'colormap')
  cfg.colormap = getcolormap(numvals, cfg);
end
if ~isequal(cfg.thres,[0 0]) % THRESHOLD a given map
  cfg.colormap = threscolormap(cfg);
end

%% -- Figure color & base color
if ~isfield(cfg,'figurecolor') && ~isfield(cfg,'basecolor')
  if isequal(cfg.thres, [0 0])
    cfg.figurecolor = 'w';
    cfg.basecolor = 'bright';
  else
    cfg.figurecolor = 'k';
    cfg.basecolor = 'dark';
  end
end
if ~isfield(cfg,'basecolormap')
  switch cfg.basecolor
    case 'dark'
      cfg.basecolormap = [.15 .35]'*ones(1,3);
    case 'lightdark'
      cfg.basecolormap = [.5 .7]'*ones(1,3);
    case {'bright','light'}
      cfg.basecolormap = [.6 .8]'*ones(1,3);
    otherwise
      error('CFG.BASECOLOR=%s NOT RECOGNIZED!',cfg.basecolor)
  end
end

%% -- No OPENGL?
g = opengl('data');
if strcmp(g.Renderer,'None')
  if isfield(surfs,'smoothwm')
    cfg.basesurf = 'smoothwm';
  elseif isfield(surfs,'white')
    cfg.basesurf = 'white';
  end
  cfg.basecurv = false;
  % NaN values will be shown as the smallest value in the colormap:
  cfg.colormap = [cfg.basecolormap(2,:); cfg.colormap];
end

%% --- Transparent subthreshold values? or not...
if ~isfield(cfg,'subthres')
  cfg.subthres = 0;
end


%% -- Initialize figure
if ~isfield(cfg,'figurehandle')
  cfg.figurehandle = figure;
else
  clf(cfg.figurehandle)
end
set(cfg.figurehandle, 'position', cfg.figureposition);
if isfield(cfg,'fname_png') % if fname_png is given, make it invisible
  set(cfg.figurehandle, 'visible','off')
end


%% M A I N ================================================================
%% --- Plot surfaces
H = struct();
for iaxes = 1:numel(cfg.hemis)
  hemi = cfg.hemis(iaxes);
  axespos(cfg.surfaxes, iaxes);
  cfg.view = cfg.views{iaxes};
  % to show partial surfaces like supratemporal planes:
  idx = strfind(cfg.basesurf,'.');
  if ~isempty(idx)
    surf2show = surfs.(cfg.basesurf(1:idx-1)).(cfg.basesurf(idx+1:end)){hemi};
    data2show = data{hemi}(surf2show.vert_idx1_orig,:);
    cfg.mask = cfg.masks{hemi}(surf2show.vert_idx1_orig);
    cfg.curv = surfs.curv{hemi}(surf2show.vert_idx1_orig);
  else
    surf2show = surfs.(cfg.basesurf){hemi};
    data2show = data{hemi};
    cfg.mask = cfg.masks{hemi};
    cfg.curv = surfs.curv{hemi};
  end
  if isfield(cfg, 'isocluscolors') % a unique color table for each hemi
    cfg.isocluscolor = cfg.isocluscolors{hemi};
  end
  [h, cfg] = myfs_viewsurf_helper(surf2show, data2show, cfg);
  if cfg.isinteger
    h.oversurf.FaceColor = 'flat';
  end
  fn = fieldnames(h);
  for j=1:numel(fn), H(iaxes).(fn{j}) = h.(fn{j}); end
  if DEBUG
    axis on % for DEBUGGING
  end
end

%% -- Plot histograms
if isequal(get(gcf,'color'), [1 1 1])
  cfg.histlinecolor = [.5 .5 .5];
  cfg.histfontcolor = [0 0 0];
elseif isequal(get(gcf,'color'), [0 0 0])
  cfg.histlinecolor = [.5 .5 .5];
  cfg.histfontcolor = [1 1 1];
end
if ~isfield(cfg,'histfontsize')
  if ismac
    cfg.histfontsize = 9;
  elseif isunix
    cfg.histfontsize = 8;
  elseif ispc
    cfg.histfontsize = 10;
  end
end


if cfg.IsHistogram

  for hemi = 1:2
    iaxes = iaxes + 1;
    H(iaxes).axes = axes('position',cfg.histaxes{hemi});
    set(H(iaxes).axes,'color',get(gcf,'color'),'xcolor',cfg.histfontcolor, 'ycolor',cfg.histfontcolor, ...
      'fontsize',cfg.histfontsize );
    hold on;

    numvals = cat(1,data{hemi});
    numvals(~cfg.masks{hemi}) = [];
    H(iaxes).hist = colorhist(numvals, cfg);
    if isfield(cfg, 'histxtick')
      set(H(iaxes).axes, 'xtick', cfg.histxtick)
    end
    if isfield(cfg, 'histxticklabel')
      set(H(iaxes).axes, 'xticklabel', cfg.histxticklabel)
    end
    if isfield(cfg, 'histxlim')
      set(H(iaxes).axes, 'xlim', cfg.histxlim)
    end
  end
end
%% -- Plot colorbar
if cfg.IsColorbar
  numvals = cat(1,data{:});
  numvals(~cat(1,cfg.masks{:})) =[];
  numvals(isnan(numvals)) = [];
  numvals(isinf(numvals)) = [];
  if numel(unique(numvals)) ~= 2
    numvals(numvals==0) = [];
  end
  if ~isempty(numvals)
    iaxes = iaxes + 1;
    H(iaxes).axes = axes('position',cfg.colorbaraxes);
    if ~isfield(cfg,'colorbarfontsize')
      if ismac
        cfg.colorbarfontsize = 13;
      elseif isunix
        cfg.colorbarfontsize = 10;
      elseif ispc
        cfg.colorbarfontsize = 11;
      end
    end
    H(iaxes).colorbar = colorbar(H(iaxes).axes, 'location','north', ...
      'fontsize',cfg.colorbarfontsize*0.9);
    clim(H(iaxes).axes, cfg.caxis);
    if cfg.isinteger && (numel(unique(numvals))<20)
      % This is to descritize the colorbar for integer values:
      numbins = min([256 numel(unique(numvals))]);
      if exist('histcounts','file')
        [~, edges] = histcounts(numvals, numbins, 'BinMethod','integers');
        xi = 0.5 * (edges(1:end-1) + edges(2:end)); % center value of bins
      else
        [~, xi] = hist(numvals, numbins);
      end
      ind = rescale(xi, 0, (size(cfg.colormap,1)-1), 'InputMin',cfg.caxis(1), 'InputMax', cfg.caxis(2) );
      ind = 1+floor(ind);
      cfg.colormap = cfg.colormap(ind,:);
    end
    colormap(H(iaxes).axes, cfg.colormap)
    axis off;
    if ~isfield(cfg,'colorbarxlabelinterp')
      cfg.colorbarxlabelinterp = 'none';
    end
    if isfield(cfg,'colorbarxlabel')
      xlabel(H(iaxes).colorbar, cfg.colorbarxlabel, ...
        'fontsize', cfg.colorbarfontsize, 'color',cfg.histfontcolor, 'interpreter',cfg.colorbarxlabelinterp);
    end
    if ~isfield(cfg,'colorbarinterp')
      cfg.colorbarinterp='none';
    end
    if isfield(cfg,'colorbartitle ')
      title(H(iaxes).colorbar, cfg.colorbartitle, ...
        'fontsize', cfg.colorbarfontsize*1.2,'fontweight','bold', ...
        'color',cfg.histfontcolor, 'interp',cfg.colorbarinterp);
    end
    if ~isfield(cfg,'colorbarxtick')
      if (cfg.isinteger) && (numel(unique(numvals))<20)
        xtick = xi;
        step = (xi(end)-xi(1))/numel(xi);
        xtick2 = xtick(1) + step/2 + [0:step:step*(numel(xi)-1)];
        cfg.colorbarxtick = xtick2;
        cfg.colorbarxticklabel = xtick;
      else
        cfg.colorbarxtick = get(H(iaxes).colorbar, 'xtick');
      end
    end
    set(H(iaxes).colorbar,'color',get(gcf,'color'), ...
      'xcolor',cfg.histfontcolor, 'ycolor',cfg.histfontcolor, 'xtick',cfg.colorbarxtick);
    if isfield(cfg,'colorbarxticklabel')
      set(H(iaxes).colorbar, 'xticklabel', cfg.colorbarxticklabel)
    end
  end
end

%% -- Print
if isfield(cfg,'fname_png')
  if ~isfield(cfg,'dpi')
    cfg.dpi = 300;
  end
  if ~isempty(getCurrentWorker)
    rendopt = '-painters';
  else
    rendopt = '-opengl';
  end
  drawnow; pause(1)
  export_fig(cfg.figurehandle, cfg.fname_png, ['-r',num2str(cfg.dpi)], rendopt)
  close(cfg.figurehandle)
end
if ~nargout, clear H cfg; end

end

function TEST()

parfor i = 1:2
  surfs = myfs_readsurfs('fsaverage4');
  thn = surfs.thickness;
  thn{1}(thn{1}==0) = nan;
  thn{2}(thn{2}==0) = nan;
  myfs_viewsurf(surfs, thn, ...
    struct('fname_png',['~/test',num2str(i),'.png']))
end


end
