function compare_unwarped (fn_distorted, fn_unwarped, fn_t1w, cfg)
% compare_unwarped (fn_distorted, fn_unwarped, fn_t1w, cfg)
%
% (cc) 2019,2020, sgKIM.

if ~nargin, help(mfilename); return; end
if ~exist('cfg','var'), cfg = []; end
[p2,f2,~] = myfileparts(fn_unwarped);
fn_gif = [p2,'/',f2,'.gif'];
if ~isfield(cfg,'figureposition')
%   cfg.figureposition = get(0,'DefaultFigurePosition');
  cfg.figureposition = [1 1 150*5*3 150*5*3];
end
cfg.figurehandle = figure('color','k', 'position',cfg.figureposition,'visible','off');
cfg.contour = fn_t1w;
cfg.xyz = 'axi25';
cfg.layout = [5 5];
I = {fn_distorted, fn_unwarped};
[Img] = niftireadgz(fn_distorted);
Img = mean(double(Img),4);
cfg.basecaxis = [min(Img(:)), max(Img(:))];
Title = {'Distorted (uncorrected)','Unwarped (corrected)'};
delay_sec = [2 2];
ColorMaps = {gray, bone};
for i = 1:2
  clf
  cfg.basecolormap = ColorMaps{i};
  slices(I{i}, [], cfg)
  axes('position',[0 0 1 0.93])
  ht = title(Title{i},'color','w', 'interp','none','fontweight','normal',...
    'horizontalAlignment','left','position',[0.001 1.0092 0.5]);
  axis off
  gifani(gcf, fn_gif, i, delay_sec(i));
end
end
