function colormap = getcolormap(numvals, cfg)
% colormap = getcolormap(numvals, cfg)

if ~isfield(cfg,'islabel')
  cfg.islabel = isequal(numvals, round(numvals)) && numel(unique(numvals))<50;
end

% Set # of color levels
if cfg.islabel
  numcolorlevels = numel(unique(numvals));
else
  numcolorlevels = 256; % default # of color levels
end

% Find colormap:
if isequal(cfg.thres, [0 0])
  % UNTHRESHOLDED
  if prod(cfg.caxis)<0 % digervent
    colormap = flipud(brewermap(numcolorlevels,'spectral'));
  else
    if cfg.caxis(2) <= 0 % negative
      cmap = flipud(brewermap(numcolorlevels*2,'spectral'));
      colormap = cmap(1:numcolorlevels,:); % green to yellow
    else
      cmap = flipud(brewermap(numcolorlevels*2,'spectral'));
      colormap = cmap(numcolorlevels+1:end,:); % yellow to red
    end
  end
else
  % THRESHOLDED
  ind = linspace(cfg.caxis(1),cfg.caxis(2),numcolorlevels);
  numlevels_neg = sum(ind <= cfg.thres(1));
  numlevels_subthr = sum(cfg.thres(1) < ind & ind < cfg.thres(2));
  numlevels_pos = sum(cfg.thres(2) <= ind);
  colormap = [
    [linspace(0,0,numlevels_neg); linspace(1,0.25,numlevels_neg); linspace(1,1,numlevels_neg)]'
    0.35*ones(numlevels_subthr,3)
    [linspace(1,1,numlevels_pos); linspace(0.25,1,numlevels_pos); linspace(0,0,numlevels_pos)]'
    ];
end

end

function TEST()
getcolormap(nan,struct('thres',[0 0],'caxis',[0 1]))
end
