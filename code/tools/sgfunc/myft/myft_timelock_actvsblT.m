function [timelock_withTrials,epoched] = myft_timelock_actvsblT(CFG)
% CFG requires:
%  .fn_data  '1xN' filename for preprocessed data (.mat)
%  .prestim  [1x1] time interval before stim onset [sec]
%  .poststim [1x1] time interval after stim onset [sec]
%  .cond     (1xK) condition structure
%  .cond(k).name '1xN' condition name
%  .cond(k).keyval {1xP} key-value pairs

fn_data=CFG.fn_data;
[p1,f1,e1]=myfileparts(fn_data);
prefix=[p1,'/',f1,'_-',num2str(CFG.prestim),'-',num2str(CFG.poststim),'s'];
if ~isfield(CFG,'channel'), CFG.channel='AG*'; end

% epoch trials (quick)
[epoched,CFG] = myft_epoch(CFG);

% compute average (also quick)
timelock=[];
timelock_withTrials=[];
for k=1:numel(CFG.cond)
  name=CFG.cond(k).name;
  nTrl=numel(epoched.(name).trial);
  
  % compute average & variance
  cfg=struct('keeptrials','no', 'channel','AG*');
  timelock.(name) = ft_timelockanalysis(cfg, epoched.(name));
  
  figure;
  df=median(timelock.(name).dof(:));
  titlestr={fn_data,[name,'(df=',(num2str(df)),')']};
  myft_plot_timelockrms(timelock.(name), titlestr)
  fn_png=[prefix,'_',name,'_ERF.png'];
  screen2png(fn_png,150);
  close(gcf) 
  
  cfg=struct('keeptrials','yes', 'channel','AG*');
  ft_warning off
  timelock_withTrials.(name) = ft_timelockanalysis(cfg, epoched.(name));
end

save([prefix,'_ERF.mat'],'timelock')
end

