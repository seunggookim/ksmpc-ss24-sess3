function [surfs, lut] = fsss_isoclus(surfs, labels, cfg)
% [surfs, lut] = fsss_isoclus(surfs, labels, cots, cfg)
%
% taking `label` cell array <1x2> as input
% returning surfs.(cfg.basesurf).isoclus structure to be used in FSSS_VIEW
%
% (cc) 2020, sgKIM.

if ~exist('cfg','var')
  cfg = [];
end
if ~isfield(cfg,'basesurf')
  cfg.basesurf = 'INFL';
end
lut = cell(1,2);
for s = 1:2
  [l,~,idx] = unique(labels{s}); % convert it into continous labels from 1 to N
  if isfield(cfg,'cots')
    lut{s} = zeros(max(idx),3);
    cot = cfg.cots{s}.table;
  end
  surfs.(cfg.basesurf){s}.isoclus = struct();
  for i = 1:max(idx)
    smoothclus = SurfStatSmooth(...
      double(idx==i)', surfs.(cfg.basesurf){s}, 4); 
    % smoothing a little to beautify it a little (but now you need SurfStat)
    
    [~,c] = tricontour(surfs.(cfg.basesurf){s}, smoothclus, 0.5, false);
    surfs.(cfg.basesurf){s}.isoclus(i).group = c.group;
    
    if isfield(cfg,'cots')
      idx_l = cot(:,5) == l(i);
      if sum(idx_l)
        lut{s}(i,:) = cot(idx_l,1:3)/255;
      else
        lut{s}(i,:) = [nan nan nan]; % medial walls
      end
    end
    
  end
end