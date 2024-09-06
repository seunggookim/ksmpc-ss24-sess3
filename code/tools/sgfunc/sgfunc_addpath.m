function sgfunc_addpath()
[mypath,~,~] = fileparts(mfilename('fullpath'));
paths = genpath(mypath);
paths_cell = strsplit(paths,':');
paths_cell = paths_cell(~contains(paths_cell,'trash'));
addpath(strjoin(paths_cell,':'))
warning off MATLAB:prnRenderer:opengl
warning off export_fig:exportgraphics

% MATLAB graphic setup ----
cmap = get(0, 'FactoryAxesColorOrder');
% idx = [1 2 5 3 4 6 7];
idx = [6 2 5 3 4 1 7];
set(0,'DefaultAxesColorOrder', cmap(idx,:));
clear cmap idx
set(0, 'DefaultLegendBox','off')
set(0, 'DefaultFigureColor','w')
set(0, 'DefaultLegendLocation','best')
set(0, 'DefaultAxesFontname','Ubuntu')
set(0, 'DefaultFigurePosition',[1 300 560 420])

end
