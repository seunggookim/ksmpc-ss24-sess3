function hHist = addmarginalhist(hScatter, Cfg)
% hHist = addmarginalhist([hScatter], [Cfg])
% 2023, dr.seunggoo.kim@gmail.com

if not(exist('hScatter','var'))
  hScatter = findobj(get(gca,'Children'),'type','Scatter');
end
if not(exist('Cfg','var'))
  Cfg = [];
end
Cfg = defaultcfg(struct(HistHeight=.08, HistMargin=.1, XylineColor=.6*[1 1 1], HistLineWidth=2, ...
  HistColor=.7*[1 1 1], HistAlpha=.5), Cfg, mfilename);

X = hScatter.XData;
Y = hScatter.YData;
Xlim = xlim; Ylim = ylim;
[Fx, Ix] = ksdensity(X);
[Fy, Iy] = ksdensity(Y);

hold on
hHist = [];
nFx = Fx./max(Fx)*range(Ylim)*Cfg.HistHeight + Ylim(1) - range(Ylim)*Cfg.HistMargin; % normalize + shift
hHist(1) = patch([Ix Ix(1)], [nFx nFx(1)], 0, LineStyle='none', FaceColor=Cfg.HistColor, FaceAlpha=Cfg.HistAlpha);
nFy = Fy./max(Fy)*range(Xlim)*Cfg.HistHeight + Xlim(1) - range(Xlim)*Cfg.HistMargin;
hHist(2) = patch([nFy nFy(1)], [Iy Iy(1)], 0, LineStyle='none', FaceColor=Cfg.HistColor, FaceAlpha=Cfg.HistAlpha);
hHist(3) = xline(Xlim(1), Color=Cfg.XylineColor); 
hHist(4) = yline(Ylim(1), Color=Cfg.XylineColor);
xlim([Xlim(1)-range(Xlim)*Cfg.HistMargin, Xlim(2)])
ylim([Ylim(1)-range(Ylim)*Cfg.HistMargin, Ylim(2)])
hold off

end
