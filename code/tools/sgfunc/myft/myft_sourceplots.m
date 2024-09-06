function h = myft_sourceplots(cfg, func, anat, titletext)
% h = myft_sourceplots(cfg, func, anat, titletext)
% (cc) 2020, sgKIM.

F = double(anat.tri);
V = anat.pos;
views = [180 0; 0 90; -90 0; 0 0; 0 -90; 90 0];
viewnames = {'front','superior','left','back','inferior','right'};
ax = axeslayout([2 3],[.1, .03, .07, .13],[.005 .005 .1 .005]);
vals = func.avg.(cfg.funparameter);
clim = [-max(abs(vals(:))) max(abs(vals(:)))];

for i = 1:6
  axespos(ax, i)
  h(i) = patch('faces',F, 'vertices',V, 'facecolor','interp', ....
    'edgecolor', 'none', 'FaceVertexCData', vals, ...
    'ambientstrength',0.4, 'diffusestrength',0.8, 'specularstrength',0);
  view(views(i,:))
  title(viewnames{i})
  camlight(0,10)
  caxis(clim);
  axis image
end
if exist('titletext','var')
  bigtitle(titletext)
end
end