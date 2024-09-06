function lea_addpath()
MyPath = fileparts(mfilename('fullpath'));
addpath(MyPath)

% MATLAB UTILITIES (INTERNAL & EXTRENAL)
assert(exist('sgfunc_addpath','file'), 'SGFUNC not in MATLAB PATH!')
if ~exist('randomize_phase','file')
  sgfunc_addpath()
  warning off export_fig:exportgraphics
end

end
