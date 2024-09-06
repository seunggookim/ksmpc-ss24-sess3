function [h,H] = imagescmat(X, varargin)
FontColor = [0 0 0 .5];
h = imagesc(X, varargin{:});
H = [];
hold on
for i = 1:size(X,1)
  for j = 1:size(X,2)
    H(i,j) = text(j,i, sprintf('%.2f',X(i,j)),'fontsize',7, 'horizontalAlignment','center',...
      'verticalAlignment','middle','color',FontColor);
  end
end

% TODO: PARSE the varargin for imagesc and text?