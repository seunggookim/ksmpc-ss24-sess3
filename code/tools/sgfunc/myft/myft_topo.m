function cfg=myft_topo(cfg, data)
% cfg=myft_topo(cfg, data)
% 
% (cc) 2019, sgKIM.

if ~isfield(cfg,'layout')
  cfg.layout=ft_prepare_layout(struct('layout','yokogawa160_helmet.mat'));
%   cfg.layout = ft_prepare_layout(struct('layout','easycapM24.mat'));
end
cfg.comment='no';
if ~isfield(cfg, 'markersymbol'), cfg.markersymbol='.'; end
if ~isfield(cfg,'title'), cfg.title='off'; end
cfg.dataname=inputname(2);

if contains(data.label{1}, 'ica')
  ft_topoplotIC(cfg, data)
elseif isfield(data,'avg')
  ft_topoplotER(cfg, data)
elseif isfield(data,'stat')
  cfg.parameter='stat';
  % find masked channels
  idx_t=data.time>=cfg.xlim(1) & data.time<=cfg.xlim(2);
  idx_chan=any(data.mask(:,idx_t),2);
  cfg.highlight='on';
  cfg.highlightchannel = data.label(idx_chan);
  if ~isfield(cfg,'highlightsymbol'),  cfg.highlightsymbol = 'o'; end
  if ~isfield(cfg,'highlightsize'), cfg.highlightsize = 3; end
  ft_topoplotER(cfg, data)
else
  ft_topoplotER(cfg, data)
end

if ~isfield(cfg,'zlim')
  c0=caxis;
  caxis([-max(abs(c0)) max(abs(c0))]);
end
colormap(gca,flipud(brewermap(256,'PuOr')));
if ~nargout
  clear cfg
end
end