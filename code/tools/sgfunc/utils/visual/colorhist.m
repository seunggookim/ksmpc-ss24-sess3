function [h_hist,cfg] = colorhist(numvals, cfg)
%COLORHIST creates colored histogram
% [USAGE]
% h_hist = colorhist(numvals, cfg)
%
% [INPUTS]
% numvals [Nx1] numerica vector. non-numeric values (nan,inf) will be removed.
% cfg     (1x1) can parse:
% (.colormap)
% (.isinteger)
% (.caxis)
% (.thres)
% (.basecolormap) [2x3] (for MYFS_VIEW) used for binary curvature
% (.subthres)     [1x1] true=make subthreshold values transparent
%                       false=replace subthreshold values with
%                       cfg.basecolormap(2,:)
%
% (cc) 2019-2020, sgKIM. solleo@gmail.com

%% - C O N F I G U R A T I O N
if ~exist('cfg','var'), cfg = []; end
if ~isfield(cfg,'histlinecolor'), cfg.histlinecolor = [.5 .5 .5]; end
if ~isfield(cfg,'histfontcolor'), cfg.histfontcolor = [0 0 0]; end
if ~isfield(cfg,'histfontsize'), cfg.histfontsize = 8; end
if ~isfield(cfg,'method'), cfg.method = []; end

%% -- Color axis default:
numvals = double(numvals);
numvals(isnan(numvals)) = [];
numvals(isinf(numvals)) = [];
if isempty(numvals)
  warning('All values are masked. Skipping a histogram.')
  h_hist = [];
  return
end
if ~isfield(cfg,'isinterger')
  cfg.isinteger = isequal(numvals, round(numvals));
end
if ~isfield(cfg,'thres'), cfg.thres = [0 0]; end
if numel(cfg.thres) == 1
  cfg.thres = [-abs(cfg.thres) abs(cfg.thres)];
end
if ~isfield(cfg,'caxis')
  cfg.caxis = winsorcaxis(numvals);
end


%% -- Color map default:
numcolorlevels = 256;
if ~isfield(cfg,'subthres'), cfg.subthres = 0; end
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
  cfg.basecolormap = [.4 .4 .4; .7 .7 .7];
end


%% Create histogram
hold on
numbins = min([2^8 numel(unique(numvals))]);
if exist('histcounts','file')
  if cfg.isinteger
    [ci, edges] = histcounts(numvals, numbins, 'BinMethod','integers');
  else
    [ci, edges] = histcounts(numvals, numbins);
  end
  xi = 0.5 * (edges(1:end-1) + edges(2:end)); % center value of bins
else
  [ci, xi] = hist(numvals, numbins);
end

% color-code bars:
ind = rescale(xi, 0, (size(cfg.colormap,1)-1), ...
  'InputMin',cfg.caxis(1), 'InputMax', cfg.caxis(2) );
ind = 1+floor(ind);
rgb = cfg.colormap(ind,:) * 0.9; % a bit darker (to see bars on white bg)

% suprathrs_neg:
idx = xi <= cfg.thres(1);
if sum(idx)
  h_neg = bar(xi(idx), ci(idx), 1, ...
    'edgecolor','none', 'facecolor','flat');
  h_neg.CData = rgb(idx,:);
else
  h_neg = [];
end

% subthres:
idx = cfg.thres(1) < xi & xi < cfg.thres(2);
if sum(idx)
  h_subthres = bar(xi(idx), ci(idx), 1, ...
    'edgecolor','none', 'facecolor','flat');
  h_subthres.CData = rgb(idx,:);
else
  h_subthres = [];
end

% suprathrs_pos:
idx = cfg.thres(2) <= xi;
if sum(idx)
  h_pos = bar(xi(idx), ci(idx), 1, 'edgecolor','none', 'facecolor','flat');
  h_pos.CData = rgb(idx,:);
else
  h_pos = [];
end

try
  xlim([edges(1) edges(end)])
catch
  xlim([min(numvals) max(numvals)]);
end
ylim0 = ylim;
box on;

% color-limits:
line([cfg.caxis(1);cfg.caxis(1)],ylim0, ...
  'color', cfg.colormap(1,:)*.85, 'linewidth',0.5, 'linestyle',':');
line([cfg.caxis(2);cfg.caxis(2)],ylim0, ...
  'color', cfg.colormap(end,:)*.85, 'linewidth',0.5, 'linestyle',':');

% thresholding:
if ~isequal(cfg.thres, [0 0])
  line([cfg.thres;cfg.thres],[ylim0' ylim0'],...
    'color',cfg.histlinecolor, 'linewidth',0.5, 'linestyle',':');
  idx = cfg.thres(1) < xi & xi < cfg.thres(2);
  if any(idx)
    if cfg.subthres
      h_subthres.FaceAlpha = 0.5;
    end
  end
end
hold off

% RETURN handles of bars
h_hist = [h_neg, h_subthres, h_pos];
end
