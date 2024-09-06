function [H, cfg] = myfs_viewsurf_helper (surf, data, cfg)
% [H, cfg] = MYFS_VIEWSURF_HELPER (surf, data, cfg)
%
% visualize scaled data on a surface
%
% [INPUT]
% surf  (1x1)  MATLAB patch structure
% data  [VxP]  Vertex-mapped scalar data (vector or rgb-values; P=1 or 3)
% view  [1x2]  view angle, default = [-90 0]
%
% cfg   (1x1)  structure to configure options (optional):
% - color schemes (presets)
%  .style    '1xN'  'spectral' | 'spectraltrans' | 'fsl' | 'freesurfer'
% - color coding (manual)
%  .basecolor '1xN'  'dark' (0.15/0.35) | 'bright' (0.7/1)
%  .colormap  [Nx3]  default = [bluish; gray; redish] if .thres ~= 0;
%                    brewermap(256,'spectral') if .thres == 0
% - optional colormaps:
%  .caxis     [1x2]  default = [1%tile 99%tile] for positive values;
%                    [-/+max(99%tile of abs)] for sigend values
%  .thres     [1x1] | [1x2]  default = [0 0]
%  .subthres  [1x1]  true | false (default) - show subthreshold values
%                    transparently
%  .mask      [Vx1]  vertex-mapped logical vectors for masking
%
% - contours
% -- isocurvature contours:
%  .doisocurv        [1x1]  true (requires surf.isocurv)
%                           | false (default)
%  .isocurvcolor     [1x3]  default = [0 0 0]
%  .isocurvlinewidth [1x1]  default = 0.5
%  .isocurvlinestyle '1xN'  default = '-'
% -- cluster contours:
%  .doclusid         [1x1]  true (requires surf.isoclus) | false (default)
%  .isocluscolor     [1x3]  default = [1 1 1]
%  .isocluslinewidth [1x1]  default = 2
%  .isocluslinestyle '1xN'  default = '-'
% -- outerboundary contours:
%  .doouterboundary        [1x1]  true | false (default)
%  .outerboundarycolor     [1x3]  default = [0 0 0]
%  .outerboundarylinewidth [1x1]  default = 2
%  .outerboundarylinestyle '1xN'  default = '-'
%
% - lighting
%  .camlight    [1x1]  true (default) | false
%  .camposition [1x2]  default = [0 10]
%  .figurecolor [1x3] | '1x1'  default = 'k' if .thres ~= 0
%                                      = 'w' if .thres == 0
%
% [OUTPUT]
% H.basesurf
% H.oversurf
% H.light
% H.axes
% H.contour_isocurv
%
% (cc) 2019-2020, sgKIM. mailto://solleo@gmail.com
% This is distributed from https://github.com/solleo/surfviz
%
% SEE ALSO: FSSS_VIEW, FSSS_READ_ALL_FS_SURFS

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


%% C O N F I G ============================================================
if ~exist('cfg','var'), cfg=[]; end

%% Color axis
if islogical(data), data = uint8(data); end
numvals = data;
if isfield(cfg,'mask'), numvals(~cfg.mask) =[]; end
numvals(isnan(numvals)) = [];
numvals(isinf(numvals)) = [];

if ~isfield(cfg,'thres')
  cfg.thres = [0 0];
end
if numel(cfg.thres) == 1
  cfg.thres = [-abs(cfg.thres) abs(cfg.thres)];
end
if ~isfield(cfg,'subthres')
  cfg.subthres = false;
end
if ~isfield(cfg,'caxis')
  if numel(unique(numvals)) < 2 % if it's binary values
    cfg.caxis = [min(numvals) max(numvals)];
  else % Winsorization
    if all(numvals>=0)
      cfg.caxis = [prctile(numvals,1) prctile(numvals,99)];
    else
      cfg.caxis = [-prctile(abs(numvals),99) +prctile(abs(numvals),99)];
    end
  end
end


%% M A I N ================================================================

%% Base surface: binary curvature ([0,1]) or just [1])
if ~isfield(cfg,'curv')
  img_base = ones(size(data,1));
else
  img_base = single(cfg.curv<0);
end
if cfg.subthres
  img_base = single(cfg.curv<0);
  img_base(cfg.mask) = ones(sum(cfg.mask),1);
end
img_base_rgb = squeeze(ind2rgb(img_base+1, cfg.basecolormap));
% see CAT_SURF_RENDER for an example of using ind2rgb


%% Overlay surface
img_over = data;
img_over(cfg.thres(1)<data & data<cfg.thres(2)) = nan; % thresholding
if isfield(cfg,'mask'), img_over(~cfg.mask) = nan; end % masking
g = opengl('data');
if strcmp(g.Renderer, 'None')
  img_over(isnan(img_over)) = 0;
end


%% Surfaces
H.axes = gca;
hold on
V = surf.vertices;
F = surf.faces;
if isfield(cfg,'basecurv') && ~cfg.basecurv
  H.basesurf = [];
else
  H.basesurf = patch('faces',F, 'vertices',V, 'facecolor','interp', ...
    'edgecolor', 'none', 'FaceVertexCData', img_base_rgb, ...
    'ambientstrength',0.4, 'diffusestrength',0.8, 'specularstrength',0);
end
if isfield(cfg,'basesurfcolor')
  H.basesurf.FaceColor = cfg.basesurfcolor;
end
% see CAT_SURF_RENDER for a nice exmaple of PATCH options
if ~isfield(cfg,'facealpha'), cfg.facealpha = 1; end
H.oversurf = patch('faces',F, 'vertices',V, 'facecolor','interp', 'edgecolor', 'none', 'FaceVertexCData', img_over, ...
  'facealpha', cfg.facealpha, 'ambientstrength',0.4, 'diffusestrength',0.8, 'specularstrength',0);


%% subthreshold overlay image
if cfg.subthres
  img_over = data;
  if isfield(cfg,'mask'), img_over(~cfg.mask) = nan; end % masking
  H.subthsurf = patch('faces',F, 'vertices',V, 'facecolor','interp', 'edgecolor', 'none', ...
    'FaceVertexCData', img_over, 'facealpha',0.5, 'ambientstrength',0.4, 'diffusestrength',0.8, 'specularstrength',0);
end
axis off, axis tight, axis image
hold off
if isfield(cfg, 'facecolor')
  set(H.oversurf, 'FaceColor',cfg.facecolor);
end
if ~isfield(cfg, 'view'),  cfg.view = [-90 0]; end
view(cfg.view)
set(gca, 'colormap',cfg.colormap)
caxis(gca, cfg.caxis);
set(gcf,'color',cfg.figurecolor)


%% Lighting
if ~isfield(cfg,'camlight'), cfg.camlight = 1; end
if ~isfield(cfg,'camposition'), cfg.camposition=[0 10]; end
if cfg.camlight
  if ismac && isempty(getCurrentWorker) % only MAC but not in parallel
    %===lines from CAT_SURF_RENDER (CAT12/SPM toolbox)===
    % set inner light
    H.light = light('Position',[0 0 0]);
    set(H.basesurf,'BackFaceLighting','unlit');
    set(H.oversurf,'BackFaceLighting','unlit');
    if isfield(H,'subthsurf')
      set(H.subthsurf,'BackFaceLighting','unlit');
    end
    %===lines from CAT_SURF_RENDER (CAT12/SPM toolbox)===
  else
    H.light = camlight(cfg.camposition(1), cfg.camposition(2));
  end
end


%% Isocurvature contours
if ~isfield(cfg,'isocurvcolor')
  cfg.isocurvcolor = [0 0 0];
end
if ~isfield(cfg,'isocurvlinewidth')
  cfg.isocurvlinewidth = 0.5;
end
if ~isfield(cfg,'isocurvlinestyle')
  cfg.isocurvlinestyle = '-';
end
if isfield(surf,'isocurv') && ~(isfield(cfg,'doisocurv') && ~cfg.doisocurv)
  hold on
  H.contour_isocurv = hggroup;   % to group LOTS of lines
  for g = 1:numel(surf.isocurv)
    xyz = surf.isocurv(g).xyz;
    plot3(xyz(:,1), xyz(:,2), xyz(:,3), 'color',cfg.isocurvcolor, ...
      'linestyle',cfg.isocurvlinestyle,'linewidth',cfg.isocurvlinewidth, ...
      'Parent',H.contour_isocurv)
  end
  hold off
end

%% Isocluster contours
if isfield(surf,'isoclus')
  if ~isfield(cfg,'isocluscolor')
    cfg.isocluscolor = repmat([1 1 1], [numel(surf.isoclus),1]);
  end
  if numel(surf.isoclus) ~= size(cfg.isocluscolor,1)
    cfg.isocluscolor = repmat(cfg.isocluscolor(1,:), [numel(surf.isoclus),1]);
  end
end
if ~isfield(cfg,'isocluslinewidth')
  cfg.isocluslinewidth = 2;
end
if ~isfield(cfg,'isocluslinestyle')
  cfg.isocluslinestyle = '-';
end
if isfield(surf,'isoclus') && ~(isfield(cfg,'doisoclus') && ~cfg.doisoclus)
  hold on
  H.contour_isoclus = hggroup;  % to group LOTS of lines
  for c = 1:numel(surf.isoclus)
    for g = 1:numel(surf.isoclus(c).group)
      xyz = surf.isoclus(c).group(g).xyz;
      if ~isnan(cfg.isocluscolor(c,:))
      plot3(xyz(:,1), xyz(:,2), xyz(:,3), 'color',cfg.isocluscolor(c,:), ...
        'linestyle',cfg.isocluslinestyle,'linewidth',cfg.isocluslinewidth, ...
        'Parent',H.contour_isoclus)
      end
    end
  end
  hold off
end


%% Outerboundary contours
if ~isfield(cfg,'outerboundarycolor')
  cfg.outerboundarycolor = [0 0 0];
end
if ~isfield(cfg,'outerboundarylinewidth')
  cfg.outerboundarylinewidth = 0.5;
end
if ~isfield(cfg,'outerboundarylinestyle')
  cfg.outerboundarylinestyle = '-';
end
if isfield(surf,'outerboundary') && ~(isfield(cfg,'doouterboundary') ...
    && ~cfg.doouterboundary)
  hold on
  H.contour_outerboundary = hggroup;
  for g = 1:numel(surf.outerboundary)
    xyz = surf.outerboundary.group(g).xyz;
    plot3(xyz(:,1), xyz(:,2), xyz(:,3), ...
      'color',cfg.outerboundarycolor, ...
      'linestyle',cfg.outerboundarylinestyle, ...
      'linewidth',cfg.outerboundarylinewidth, ...
      'Parent',H.contour_outerboundary)
  end
  hold off
end

end
