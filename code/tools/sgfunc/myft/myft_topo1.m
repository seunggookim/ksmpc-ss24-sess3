function myft_topo1(y, refData, cfg)
% quick view of a given vector (y) and refData structure;
%
% myft_topo1(y, refData, cfg)
%
dat =[];
if exist(dat,'grad')
  dat.grad = refData.grad;
end
dat.label = refData.label;
dat.dimord = 'chan_time';
dat.time = 1;
dat.y = reshape(y,[],1);
if ~exist('cfg','var')
  cfg=[];
end
cfg.xlim = [1 1];
cfg.parameter = 'y';
% cfg.contours = 0;
% cfg.markersymbol = 'o';
% cfg.contournum = 0;
myft_topo(cfg, dat);