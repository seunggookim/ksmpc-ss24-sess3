function myspm_plot(cfg)
% myspm_plot(cfg)

load (cfg.fname_plotdata);
figpos=[1924         564        126.8929+129.6429         201];
figure('color','k', 'position',figpos);
axes('position',[0.25 0.24 .65 .6]);
hold on;
line(plotdata.x, plotdata.yhat,'color',[1 0 0],'linewidth',2)
scatter(plotdata.x, plotdata.y,'wo')
xlabel(cfg.xlabel,'color','w','fontsize',11);
ylabel(cfg.ylabel,'color','w','fontsize',11);
set(gca,'color','k','xcolor','w','ycolor','w')
% title label
if ~isfield(cfg,'titlefontsize'), cfg.titlefontsize=13; end
if isfield(cfg,'title')
 hal=axes('position',[0 0.94 1 0.1]); axis off;
 if isfield(cfg,'titlepositiony')
  set(hal,'position',[0 cfg.titlepositiony 1 0.1]);
 end
 text(hal,0,0, cfg.title, 'fontsize',cfg.titlefontsize, ...
 'color','w', 'interp','none', 'BackgroundColor', 'none', ...
 'HorizontalAlignment', 'left', 'VerticalAlignment','Middle');
end

if ~isfield(cfg,'dpi')
 dpi=900;
else
 dpi=cfg.dpi;
end
if isfield(cfg,'fname_png')
 screen2png(cfg.fname_png,dpi);
end
end