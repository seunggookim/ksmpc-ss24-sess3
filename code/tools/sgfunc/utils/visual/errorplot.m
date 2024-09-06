function [h,means,sem]=errorplot(x,Y,se,useErrorBar)
% h=errorplot(x,Y,se)
% draws plots with SE (default=1)
%
% y is a 2-D matrix in [# of samples] x [# of conditions (or # of x-bins)]
%
% (cc) 2015, sgKIM. solleo@gmail.com
if nargin==0, help(mfilename); return; end
if ~exist('se','var')
  se = 1;
end
if isempty(x)
  x = 1:size(Y,2);
end
if size(Y,1)~=1 && size(Y,2)==1
  Y = Y';
end
means = mean(Y,1);
sem = se * std(Y,[],1) ./ sqrt(size(Y,1));
if ~exist('useErrorBar','var'), useErrorBar=false; end
if exist('shadedErrorBar.m','file') && ~useErrorBar
  h = shadedErrorBar(x,means,sem);
else
  h = errorbar(x,means,sem);
end
if ~nargout, clear h means sem; end
end
