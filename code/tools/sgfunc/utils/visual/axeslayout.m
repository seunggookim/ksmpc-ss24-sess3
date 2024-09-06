function ax = axeslayout (layout,inset,outset)
% ax = axeslayout (layout,inset,outset)
%
% layout = [RxC] # of rows, # of columns of a grid of plotboxes to construct
% inset  = [1x4] <left right top bottom> inset margin within the grid
%                 values between 0 and 1 w.r.t a plotbox size
%          or 'tight' [.1, .03, .07, .13]
% outset = [1x4] <LEFT RIGHT TOP BOTTOM> outset margin outside the grid
%                 values between 0 and 1 w.r.t a figure size
%          or 'tight' for [.005 .005 .005 .005]
%
% -------------------------------------------------
%                        TOP
% -------------------------------------------------
%      |....................................|
%      |:      top       ::      top       :|
%      |....................................|
% LEFT |:left:(^o^):right::left:(TxT):right:| RIGHT
%      |....................................|
%      |:     bottom     ::     bottom     :|
%      |....................................|
% -------------------------------------------------
%                      BOTTOM
% -------------------------------------------------
%
% Example:
% ax = axeslayout([2,4],'tight','tight');  % creates 2 rows x 4 columns axes
% figures
% for i = 1:8
%   axespos(ax, i)
%   title(sprintf('Axes %i',i))
% end
%
% SEE ALSO: AXESPOS
% (cc) sgKIM, 2019.

if ~exist('inset','var') || isempty(inset)
  inset = [.25, .03, .14, .25];  % [left right top bottom]
elseif strcmp(inset,'tight')
  inset = [.1, .03, .07, .13];  % [left right top bottom]
end
if ~exist('outset','var') || isempty(outset)
  outset = [.02, .02, .02, .02]; % [LEFT RIGHT TOP BOTTOM]
elseif strcmp(outset,'tight')
  outset = [.005 .005 .005 .005]; % [LEFT RIGHT TOP BOTTOM]
end

%% contruct a grid
W = (1-sum(outset(1:2)))/layout(2);
H = (1-sum(outset(3:4)))/layout(1);
if layout(2) > 1
  X = repmat(outset(1):W:(1-W-outset(2)), [1,layout(1)]);
else
  X = repmat(outset(1), [1,layout(1)]);
end
if layout(1) > 1
  Y = reshape(repmat((1-H-outset(3)):-H:outset(4), [layout(2),1]),[],1);
else
  Y = repmat(outset(4), [1,layout(2)]);
end

debug=0;
if debug
  plotbox=struct('x',X,'y',Y,'w',W,'h',H);
  clf
  for k=1:prod(layout)
    axespos(plotbox,k)
    set(gca,'xtick',[],'ytick',[])
  end
end

w = W * (1-sum(inset(1:2)));
h = H * (1-sum(inset(3:4)));
x = X + W*inset(1);
y = Y + H*inset(4);
ax = struct('x',x,'y',y,'w',w,'h',h);
if debug
  hold on
  for k=1:prod(layout)
    axespos(axes_pos,k)
    plot(rand(10)*10000)
    xlabel({'OOOO'}); ylabel({'YYYY'}); legend({'XXXXXXXX'},'location','east')
    title({'XXXXXXXX'})
  end
end

end