function [Z, XbinEdges, YbinEdges, normZ] = jhist(X, Y, N, XbinEdges, YbinEdges)
% [Z, Xbin, Ybin, normZ] =jhist(X, Y, N, Xbin, Ybin)
% N = # of bins (default=256)
% sgKIM, 2013-09-20.

if ~exist('N','var')
  N=[256 256];
else
  if numel(N)==1
    N=[N N];
  end
end
if nargin<=3
  XbinEdges=linspace(min(X),max(X),N(1)+1);
  YbinEdges=linspace(min(Y),max(Y),N(2)+1);
else
  N=[numel(YbinEdges) numel(XbinEdges)];
end
Z=zeros(N(2),N(1));
for x=1:numel(XbinEdges)-1
  if x==1
    ind= (XbinEdges(x) <= X) & (X <= (XbinEdges(x+1)));
  else
    ind= (XbinEdges(x) < X) & (X <= (XbinEdges(x+1)));
  end
  Z(:,x)=histcounts(Y(ind), YbinEdges); % row as Y-axis, column as X-axis
end
normZ=Z./numel(X);
if nargout == 0
  imagesc(XbinEdges, YbinEdges, Z)
%     contourf(Z)
  set(gca,'ydir','nor')
  clear Z Xbin Ybin
end
end
