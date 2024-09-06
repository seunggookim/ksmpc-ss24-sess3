function [H, cfg] = orthoslices(base, data, cfg)
if ~exist('data', 'var'), data = []; end
if ~exist('cfg', 'var'), cfg = struct(); end
cfg_slices = cfg;
if isfield(cfg, 'fname_png') % if fname_png is given, make it invisible
  cfg_slices.figurehandle = figure('visible','off');
  cfg_slices = rmfield(cfg_slices, 'fname_png');
end
if ~isfield(cfg,'xyz')
  cfg_slices.xyz = ['sag1';'cor1';'axi1'];
end
if ischar(data)
  isFsl = contains(data, 'fsl');
else
  isFsl = false;
end
[H, cfg_slices, base, data] = slices(base, data, cfg_slices);
cfg.figurehandle = cfg_slices.figurehandle;
H(3).colorbar.Position(2) = 0.475; % move the colorbar a bit above

%% Create the 4th axes
H(4) = H(3);
Fieldnames = fieldnames(H(4));
for i = 1:numel(Fieldnames)
  H(4).(Fieldnames{i}) = [];
end

pos = [H(2).baseaxes.Position(1) H(3).baseaxes.Position(2) ...
  H(2).baseaxes.Position(3:4)];
MarginWidth = pos(3)*0.15;
MarginHeight = pos(4)*0.15;
H(4).baseaxes = axes('position', ...
  [pos(1)+2*MarginWidth pos(2)+2*MarginHeight ...
  pos(3)-3*MarginWidth pos(4)-3*MarginHeight], ...
  'color', cfg_slices.figurecolor);

%% Take numerical values of...
if isempty(data)
  % BASE
  numvals = base.vol(:);
  numvals(isnan(numvals)) = [];
%   caxis = cfg_slices.basecaxis;
elseif isFsl
  ijk0_contour = find3(cfg.contour.vol>0)-1;
  ijk0_contour = [ijk0_contour ones(size(ijk0_contour,1),1)];
  xyz_contour = ijk0_contour*cfg.contour.vox2ras';
  
  ijk0_data = round(xyz_contour/data.vox2ras');
  ijk1_data = ijk0_data(:,1:3)+1;
  ijk1_grid = find3(ones(size(data.vol)));
  idx = ismember(ijk1_grid, ijk1_data,'rows');
  
  numvals = data.vol(idx);
  numvals(isnan(numvals)) = [];
  numvals(isinf(numvals)) = [];
%   caxis = cfg_slices.caxis;
else
  % DATA
  numvals = data.vol(:);
  numvals(~cfg_slices.mask(:)) = [];
  numvals(isnan(numvals)) = [];
  numvals(isinf(numvals)) = [];
%   caxis = cfg_slices.caxis;
end

%% CREATE a histogram
% H(4).baseslice = histogram(numvals, 'BinLimits', caxis);
% H(4).baseslice.FaceColor = 'c';
% H(4).baseslice.EdgeColor = 'none';
H(4).baseslice = colorhist(numvals, cfg_slices);

set(H(4).baseaxes, 'color', cfg_slices.figurecolor)
H(4).baseaxes.XColor = cfg_slices.coordfontcolor;
H(4).baseaxes.YColor = cfg_slices.coordfontcolor;

if isfield(cfg_slices, 'histtitle')
  title(H(4).baseaxes, cfg_slices.histtitle, ...
    'color',cfg_slices.coordfontcolor, 'interpret','none')
end
if isfield(cfg_slices, 'histxlabel')
  xlabel(H(4).baseaxes, cfg_slices.histxlabel, ...
    'color',cfg_slices.coordfontcolor, 'interpret','none')
end
if isfield(cfg_slices, 'histylabel')
  ylabel(H(4).baseaxes, cfg_slices.histylabel, ...
    'color',cfg_slices.coordfontcolor, 'interpret','none')
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
