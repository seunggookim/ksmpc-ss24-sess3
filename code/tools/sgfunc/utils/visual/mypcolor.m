function h=mypcolor(x,y,C, varargin)
% MYPCOLOR(xAxisVector,yAxisVector,MatrixToShow, cfg)
%
% (cc) 2019, sgKIM
% if ~exist('vargin','var'), varargin={}; end
% hold on
if ~nargin, help (mfilename); end
% if ~exist('cfg','var'), cfg = []; end
h = surface(x,y,-1+C*0,C,'linestyle','none',varargin{:});
% if isfield(cfg,'yscale'); set(gca,'yscale',cfg.yscale); end
% axis1 = [min(x) max(x) min(y) max(y)];
% if min(x)==max(x)
%   axis1(1:2) = [min(x)-0.5 min(x)+0.5];
% end
% if min(y)==max(y)
%   axis1(3:4) = [min(y)-0.5 min(y)+0.5];
% end
% axis(axis1)
% axis0 = axis;
% if isfield(cfg,'xtick'); set(gca,'xtick',cfg.xtick); end
% if isfield(cfg,'ytick'); set(gca,'ytick',cfg.ytick); end
% xtick0 = get(gca,'xtick');
% ytick0 = get(gca,'ytick');
% xtick0 = setdiff(xtick0,axis0([1 2]));
% ytick0 = setdiff(ytick0,axis0([3 4]));
% hold on
% if ~(isfield(cfg,'grid') && cfg.grid==0)
%   gridColor = [0.5 0.5 0.5];
%   line([xtick0; xtick0],repmat(axis0([3 4])',[1 numel(xtick0)]), ...
%     'linewidth',0.1,'color',gridColor, 'linestyle',':')
%   line(repmat(axis0([1 2])',[1 numel(ytick0)]),[ytick0; ytick0], ...
%     'linewidth',0.1,'color',gridColor, 'linestyle',':')
% end
set(gca,'ydir','nor','box','on','boxStyle','full','layer','top') % LAYER!!!
grid on
axis tight
% hold off

if ~nargout, clear h; end

end
