function [timefreq, CFG] = myft_freqanalysis(CFG)
% [timefreq, CFG] = myft_freqanalysis(CFG) - ft_freqanalysis wrapper
% CFG requires:
%  .fn_data  '1xN' filename for preprocessed data (.mat)
%  .prestim  [1x1] time interval before stim onset [sec]
%  .poststim [1x1] time interval after stim onset [sec]
%  .cond     (1xK) condition structure
%  .cond(k).name '1xN' condition name
%  .cond(k).keyval {1xP} key-value pairs
% (.prefix)
% (.remove_dss1)
%
%  .method
% (.taper)
%  .foi
% (.t_ftimwin)
% (.dt
% (.toi)
% (.pad)
% (.savemat)
%
% (cc) sgKIM, 2019.

fn_data = CFG.fn_data;
[p1,f1,e1] = myfileparts(fn_data);
prefix = [p1,'/',f1,'_-',num2str(CFG.prestim),'-',num2str(CFG.poststim),'s'];
if isfield(CFG,'zthres_dvnt') && isfield(CFG,'zthres_jmpy')
  prefix = [prefix,'_',num2str(CFG.zthres_dvnt),'_',num2str(CFG.zthres_jmpy)];
end
if isfield(CFG,'desc')
  prefix = [prefix,'_',CFG.desc];
end
CFG.prefix = prefix;
if ~isfield(CFG,'channel'), CFG.channel='AG*'; end
if ~isfield(CFG,'figures'), CFG.figures=1; end
if ~isfield(CFG,'overwrite'), CFG.overwrite=0; end
if ~isfield(CFG,'basecorrect'), CFG.basecorrect=1; end
if ~isfield(CFG,'savemat'), CFG.savemat=0; end
conds = {CFG.cond.name};
nConds = numel(conds);

% epoch trials + timelock proc
CFG2 = CFG;
CFG2.figures = 0;
CFG2.savemat = 0;
[timelock, ~, epoched] = myft_timelock(CFG2);

% compute TFR (30s)
if ~isfield(CFG,'pad')
  CFG.pad = 'nextpow2' ;        % padding method
end
tic
timefreq = [];
for k = 1:nConds
  if CFG.savemat
    fn_mat = [CFG.prefix,'_TFR_',conds{k},'.mat'];
    if ~exist(fn_mat,'file')
      tfr = ft_freqanalysis(CFG, epoched.(conds{k}));
      save(fn_mat,'tfr')
    end
    timefreq.(conds{k}) = fn_mat;
  else
    ft_warning off;
    timefreq.(conds{k}) = ft_freqanalysis(CFG, epoched.(conds{k}));
  end
end

%% plotting
if CFG.figures
  %%
  maxGlobalPower = 0; maxItpc = 0;
  if isfield(CFG,'savefig') && CFG.savefig
    figvis=0;
  else
    figvis=1;
  end
  figure('position',[1 608 150*nConds 240], 'visible',figvis)
  ax = axeslayout([2,nConds],[.03 .03 .14 .25],[0.02/12*nConds, .03/12*nConds, .02, .02]);
  h_ax = zeros(1,nConds);
  for k=1:nConds
    if ischar(timefreq.(conds{k}))
      load(timefreq.(conds{k}),'tfr')
    else
      tfr = timefreq.(conds{k});
    end
    F = tfr.fourierspctrm;
    h_ax(k) = axespos(ax,k);
    phi = abs(F);
    globalPower = squeeze(mean(mean(phi,1),2));
    cfgP=struct('yscale','log');
    mypcolor(tfr.time, tfr.freq, globalPower, cfgP);
    caxis([0 maxGlobalPower]);
    title(conds{k});
    xlabel('Time [s]'); 
    if k==1 
      ylabel('Freq [Hz]')
    end
    %     colormap(gca, sgcolormap('HSV-soft',8))
    if k==nConds
      hc=colorbaraxes('location','eastOutside');
      ylabel(hc,'Global power');
    end
    
    h_ax(k+nConds) = axespos(ax,k+nConds);
    itpc = squeeze( abs( sum(F./abs(F),1) ) / size(F,1) );
    itpc = squeeze( mean(itpc) ); % average over channels
    mypcolor(tfr.time, tfr.freq, itpc, cfgP);
    xlabel('Time [s]'); 
    if k==1
      ylabel('Freq [Hz]')
    end
    if k==nConds
      hc=colorbaraxes('location','eastOutside');
      ylabel(hc,'Global ITPC');
    end
    
    maxGlobalPower = max([maxGlobalPower max(globalPower)]);
    maxItpc = max([maxItpc max(itpc)]);
  end
  
  for k=1:nConds
    caxis(h_ax(k), [0 maxGlobalPower])
    caxis(h_ax(k+nConds), [0 maxItpc])
  end
  colormap parula
  if figvis==0
    screen2png([prefix,'_TFR.png'],150)
    close(gcf)
  end

end

end

