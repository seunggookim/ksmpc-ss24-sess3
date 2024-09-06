function [timelock_withTrials,CFG,epoched,timelock] = myft_timelock_eeg24(CFG)
% [timelock_withTrials,CFG,epoched,timelock] = myft_timelock(CFG) - ft_timelock wrapper
%
% CFG requires:
%  .fn_data  '1xN' filename for preprocessed data (.mat)
%  .prestim  [1x1] time interval before stim onset [sec]
%  .poststim [1x1] time interval after stim onset [sec]
%  .cond     (1xK) condition structure
%  .cond(k).name '1xN' condition name
%  .cond(k).keyval {1xP} key-value pairs
%
% (cc) sgKIM, 2019.

fn_data=CFG.fn_data;
[p1,f1,e1]=myfileparts(fn_data);
prefix=[p1,'/',f1,'_-',num2str(CFG.prestim),'-',num2str(CFG.poststim),'s'];
if isfield(CFG,'zthres_dvnt') && isfield(CFG,'zthres_jmpy')
  prefix=[prefix,'_',num2str(CFG.zthres_dvnt),'_',num2str(CFG.zthres_jmpy)];
end
if isfield(CFG,'desc')
  prefix = [prefix,'_',CFG.desc];
end
CFG.prefix = prefix;
if ~isfield(CFG,'prefix'), CFG.prefix = prefix; end
if ~isfield(CFG,'channel'), CFG.channel='all'; end
if ~isfield(CFG,'figures'), CFG.figures=1; end
if ~isfield(CFG,'basecorrect'), CFG.basecorrect=1; end
if ~isfield(CFG,'savemat'), CFG.savemat=1; end

% epoch trials (quick)
[epoched,CFG] = myft_epoch(CFG);

% separate the epoch into baseline and activation periods:
if ~isfield(CFG,'actvsbsl'), CFG.actvsbsl=0; end
if CFG.actvsbsl
  for k=1:numel(CFG.cond)
    actname=CFG.cond(k).name;
    bslname=[CFG.cond(k).name,'bsl'];
    epoched.(bslname)=epoched.(actname);
    cfg=struct('toilim',[0 CFG.poststim]);
    epoched.(actname) = ft_redefinetrial(cfg, epoched.(actname));
    cfg=struct('toilim',[-CFG.prestim 0 ]);
    epoched.(bslname) = ft_redefinetrial(cfg, epoched.(bslname));
    epoched.(bslname).time = epoched.(actname).time;
  end
end

% compute average (also quick)
timelock=[];
timelock_withTrials=[];
names=fieldnames(epoched);
for k=1:numel(names)
  name=names{k};
  
  % baseline correction of each trials before computing average & variance
  if CFG.basecorrect
    nTrl=numel(epoched.(name).trial);
    for t=1:nTrl
      begsample=find(epoched.(name).time{t}>=-CFG.prestim,1,'first');
      endsample=find(epoched.(name).time{t}<=0,1,'last');
      epoched.(name).trial{t} = ft_preproc_baselinecorrect(...
        epoched.(name).trial{t}, begsample, endsample);
    end
  end
  
  % compute average & variance
  cfg=struct('keeptrials','no', 'channel','all');
  timelock.(name) = ft_timelockanalysis(cfg, epoched.(name));
  
  if CFG.figures
    figure;
    df=median(timelock.(name).dof(:));
    titlestr={fn_data,[name,' (df=',(num2str(df)),')']};
    myft_plot_timelockrms(timelock.(name), titlestr)
    fn_png=[prefix,'_',name,'_ERP.png'];
    screen2png(fn_png,80,0);
    pause(1)
    close(gcf)
  end
  
  cfg = struct('keeptrials','yes', 'channel','all');
  ft_warning off
  timelock_withTrials.(name) = ft_timelockanalysis(cfg, epoched.(name));
end
if CFG.savemat
  save([prefix,'_ERP.mat'],'timelock')
end
end

