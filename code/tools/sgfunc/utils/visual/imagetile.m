function imagetile(Fnames, FnameOut, cfg)
%IMAGETILE creates a tile of images
% imagetile(Fnames, FnameOut, cfg)
%
% (cc) 2022, dr.seunggoo.kim@gmail.com

if ~exist('cfg','var'), cfg=[]; end
N = numel(Fnames);
if ~isfield(cfg,'Layout')
  cfg.Layout = [1 1]*ceil(sqrt(N));
end

FnRows = {};
for iRow = 1:cfg.Layout(1)
  idx = cfg.Layout(2)*(iRow-1)+[1:cfg.Layout(2)];
  idx(idx>N) = [];
  if numel(idx)
  FnRows{iRow} = [tempname,'.png'];
  imageconcat(Fnames(idx), FnRows{iRow}, 2, cfg);
  end
end
cfg.isPaddingLeft = true;
cfg.bgval = 0;
imageconcat(FnRows, FnameOut, 1, cfg)
end
