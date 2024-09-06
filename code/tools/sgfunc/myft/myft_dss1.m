function timelock = myft_dss1(CFG, timelock)
% timelock = myft_dss1(CFG, timelock)
%
% CFG requires:
%  .fn_data  '1xN' filename for preprocessed data (.mat)
%  .prestim  [1x1] time interval before stim onset [sec]
%  .poststim [1x1] time interval after stim onset [sec]
%  .cond     (1xK) condition structure
%  .cond(k).name '1xN' condition name
%  .cond(k).keyval {1xP} key-value pairs
% (.desc)
% (.idx_intrasubjecttemplate)
% (.intersubjecttemplate)
% (.figures)
%
% (cc) 2019, sgKIM.

fn_data=CFG.fn_data;
[p1,f1,e1]=myfileparts(fn_data);
prefix=[p1,'/',f1,'_-',num2str(CFG.prestim),'-',num2str(CFG.poststim),'s'];
if isfield(CFG,'zthres_dvnt') && isfield(CFG,'zthres_jmpy')
  prefix=[prefix,'_',num2str(CFG.zthres_dvnt),'_',num2str(CFG.zthres_jmpy)];
end
if ~isfield(CFG,'savemat'), CFG.savemat=1; end
if isfield(CFG,'desc')
  prefix=[prefix,'_',CFG.desc];
end
names=fieldnames(timelock);
K=numel(names);
for k=1:K
  if isfield(CFG,'toilim')
    timelock.(names{k}) = ft_redefinetrial(CFG, timelock.(names{k}));
  end
  time_chan_trl=permute(timelock.(names{k}).trial,[3 2 1]);
  [num_time, num_chan, num_trl] = size(time_chan_trl);
  [todss,fromdss,ratio,pwr]=dss1(time_chan_trl); % keep all components
  w=[];
  for l=1:num_trl
    w(:,l) = time_chan_trl(:,:,l)*todss(:,1);
  end
  
  timelock.(names{k}) = rmfield(timelock.(names{k}),'trial');
  timelock.(names{k}).dss_w  = mean(w,2);
  timelock.(names{k}).dss_s  = todss(:,1);
  timelock.(names{k}).dss_ratio = ratio;
  
  timelock.(names{k}).todss = todss;
  timelock.(names{k}).fromdss = fromdss;
  
  % trials filtered by DSS1
  timelock.(names{k}).dimord = 'rpt_time';
  timelock.(names{k}).trial = zeros(num_trl, num_time);
  for l=1:num_trl
    timelock.(names{k}).trial(l,:) = ...
      squeeze(time_chan_trl(:,:,l)) * todss(:,1);
  end
end
if CFG.savemat
  save([prefix,'_ERF_dss1.mat'],'timelock');
end

if isfield(CFG,'figures') && CFG.figures
  ax=axeslayout([3 K]);
  figure('position',[-1919         410        1906         542])
  lay=ft_prepare_layout(struct('layout','yokogawa160_helmet.mat'));
  ft_warning off
  for k=1:K
    tl=timelock.(names{k});
    axespos(ax,k)
    stem(tl.dss_ratio,'.')
    xlabel('Order');
    if k==1, ylabel('Power ratio'); end
    title(names{k})
    
    axespos(ax,k+K)
    plot(tl.time, tl.dss_w)
    xlabel('Time [s]');
    if k==1, ylabel('Temporal w'); end
    
    axespos(ax,k+K*2)
    ft_warning off
    tl.dimord = 'chan_time';
    myft_topo(struct('parameter','dss_s','xlim',[tl.time(1) tl.time(1)],...
      'layout',lay),tl)
  end
  screen2png([prefix,'_ERF_dss1.png'],100,1)
  pause(1)
  close(gcf)
end

end