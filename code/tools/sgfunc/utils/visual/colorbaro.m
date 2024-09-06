function cb = colorbaro(varargin)
ax=get(gca,'position');
cb=colorbar(varargin{:});
%pos=get(cb,'position');
set(gca,'position',ax);
if ~nargout
  clear cb
end
end