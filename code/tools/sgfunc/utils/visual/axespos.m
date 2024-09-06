function h = axespos(ax,i)
% axes_handle = axespos(ax, i)
%
% ax = <1xN Struct> output from AXESLAYOUT
% i  = [1x1] 1-based integer to make the i-th axes of the grid.
%
% Example:
% ax = axeslayout([2,4],'tight','tight');  % creates 2 rows x 4 columns axes
% figures
% for i = 1:8
%   axespos(ax, i)
%   title(sprintf('Axes %i',i))
% end
%
% SEE ALSO: AXESLAYOUT
% (cc) sgKIM, 2019.

if numel(ax.w)==1
  h = axes('position',[ax.x(i), ax.y(i), ax.w, ax.h]);
else
  h = axes('position',[ax.x(i), ax.y(i), ax.w(i), ax.h(i)]);
end
if nargout == 0
  clear h
end
end