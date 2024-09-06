function myft_itc(CFG)
% myft_itc(CFG) - a wrapper for my_ft_compute_itc
%
% see also MY_FT_COMPUTE_ITC
% (cc) 2019, sgKIM

itc_w = []; 
itc_x = [];
for iBand = 1:length(CFG.bpfBands)
  % epoching with filtering
  CFG.bpfreq = CFG.bpfBands(iBand,:);
  k=1;
  for iLv = 1:4
    for iEx = 1:3
      CFG.cond(k).name = sprintf('lvl%dex%d',iLv,iEx);
      CFG.cond(k).keyval = {'lvl',iLv, 'ex',iEx};
      k = k + 1;
    end
  end
  [timelock, CFG] = myft_timelock(CFG);
  
  % compute ITC
  intervals = [0 1; 1 2; 2 3];
%   intervals = [0 3; 0 2; 1 3; 0 1; 1 2; 2 3; 1 3;
%      0 0.5; 0.5 1; 1 1.5; 1.5 2; 2 2.5; 2.5 3];
  for iTime = 1:length(intervals)
    [itc_w(:,:,iTime,iBand), itc_x(:,:,iTime,iBand)] ...
      = myft_compute_itc(timelock, intervals(iTime,:), CFG.itc_x_conds);
  end
end

CFG.cond=[];
k=1; iLv = 1; iEx = 1;
CFG.cond(k).name = sprintf('lvl%dex%d',iLv,iEx);
CFG.cond(k).keyval = {'lvl',iLv, 'ex',iEx};
[~, CFG, ~, timelock] = myft_timelock(CFG);

fn_mat = [CFG.prefix,'_itc.mat'];
save(fn_mat,'timelock','CFG','itc_w','itc_x')
ls(fn_mat)

end
