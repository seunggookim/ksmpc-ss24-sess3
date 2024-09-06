function h = bigtitle(titletext, titlefontcolor, titlefontsize, titleinterp)
% h = bigtitle(titletext, titlefontcolor, titlefontsize, titleinterp)
% (cc) 2019, sgKIM.

if ~exist('titlefontcolor','var')
  titlefontcolor = [0 0 0];
end
if ~exist('titlefontsize','var')
  titlefontsize = 16;
end
if ~exist('titleinterp','var')
  titleinterp = 'none';
end
axes('position',[0 0 1 1])
h = text(0.5, 0.97, titletext, 'horizontalalignment','center', ...
  'color',titlefontcolor, 'fontsize',titlefontsize, 'interp', titleinterp);
axis off
if ~nargout
  clear h
end
end