function cfg = myft_topochan(idx_chan, ref_data, cfg)
% myft_topochan - draw channel markers only
% cfg = myft_topochan(idx_chan, ref_data, cfg)
%
% (cc) 2019, sgKIM.

if ~exist('cfg','var')
  cfg = struct('style','blank', 'highlight','on', 'highlightsymbol','o',...
    'highlightsize', 3);
end
cfg.highlightchannel = ref_data.label(idx_chan);
myft_topo1(zeros(157,1), ref_data, cfg)
