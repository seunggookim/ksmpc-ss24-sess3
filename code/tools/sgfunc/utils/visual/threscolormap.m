function colormap = threscolormap(cfg)
% colormap = threscolormap(cfg)

% THRESHOLD EXISTING COLORMAP
if (prod(cfg.caxis)<0)
  % divergent:
  givenmap = cfg.colormap;
  numcolorlevels = size(givenmap,1);
  ind = linspace(cfg.caxis(1),cfg.caxis(2),numcolorlevels);
  numlevels_neg = sum(ind <= cfg.thres(1));
  numlevels_subthr = sum(cfg.thres(1) < ind & ind < cfg.thres(2));
  numlevels_pos = sum(cfg.thres(2) <= ind);
  colormap = [
    givenmap(1:numlevels_neg,:);
    0.35*ones(numlevels_subthr,3);
    givenmap(end-numlevels_pos+1:end,:)];
else % NON-DIVERGENT????
    givenmap = cfg.colormap;
  numcolorlevels = size(givenmap,1);
  ind = linspace(cfg.caxis(1),cfg.caxis(2),numcolorlevels);
  numlevels_neg = sum(ind <= cfg.thres(1));
  numlevels_subthr = sum(cfg.thres(1) < ind & ind < cfg.thres(2));
  numlevels_pos = sum(cfg.thres(2) <= ind);
  colormap = [
    givenmap(1:numlevels_neg,:);
    0.35*ones(numlevels_subthr,3);
    givenmap(end-numlevels_pos+1:end,:)];
end
end
