function [isJmpyTrl,isDvnTrl]  = myft_reject_trials(CFG, data)
% [isJmpyTrl,isDvnTrl]  = myft_reject_trials(CFG, data)

if ~isfield(CFG,'channel')
  CFG.channel='AG*'; % ONLY FOR MEG160
end

if ~isfield(CFG,'zthres_jmpy')
  zthres_jmpy=10;
else
  zthres_jmpy=CFG.zthres_jmpy;
end
cfg1=[];
cfg1.trl=CFG.trl;
zvalue=struct('cutoff',zthres_jmpy,'trlpadding',0','artpadding',0,...
  'fltpadding',0);
zvalue.channel=CFG.channel;
cfg1.artfctdef.zvalue=zvalue;
cfg1.feedback='no';
ft_warning off; warning off;
[~,art]=ft_artifact_zvalue(cfg1,data);

jmpytrl=[];
for a=1:size(art,1)
  jmpytrl=[jmpytrl; ...
    find(data.sampleinfo(:,1)<=art(a,1) & art(a,2)<=data.sampleinfo(:,2))];
end
jmpytrl=unique(jmpytrl);
isJmpyTrl=false(size(CFG.trl,1),1);
isJmpyTrl(jmpytrl)=true;

fig=0;
if fig
  [p1,f1,~]=myfileparts(CFG.fn_data);
  figure('position',[  -1863          66        1786         881])
  X=permute(cat(3,data.trial{:}),[2 1 3]); % time x chan x trl
  ax=axeslayout([8,20],[0 0 0 0],[0 0 0 0]);
  for c=1:157
    axespos(ax,c)
    hold on;
    val=(denan(squeeze(X(:,c,:)))); % along time
    plot(data.time{1}, val,'color',[0 0 0 .2])
    plot(data.time{1}, mean(val,2),'color','r','linewidth',2)
    ylim([-10 10].*1e-13)
    text(0,0,['AG',num2str(c)],'color','w')
    if c==141
      title({CFG.fn_data,'Z-over-time'}, 'color','b',...
        'horizontalalignment','left','interpreter','none');
    end
  end
  screen2jpeg([p1,'/',f1,'_raw_epoched.jpg'],72)
  pause(2)
  close(gcf)
end

if fig
  [p1,f1,~]=myfileparts(CFG.fn_data);
  figure('position',[  -1863          66        1786         881])
  X=permute(cat(3,data.trial{:}),[2 1 3]); % time x chan x trl
  ax=axeslayout([8,20],[0 0 0 0],[0 0 0 0]);
  for c=1:157
    axespos(ax,c)
    hold on;
    zval=zscore(denan(squeeze(X(:,c,:)))); % along time
    plot(data.time{1}, zval(:,~isJmpyTrl),'color','k')
    if sum(isJmpyTrl)
      plot(data.time{1}, zval(:,~~isJmpyTrl),'color','r')
    end
    ylim([-10 10])
    text(0,0,['AG',num2str(c)],'color','w')
    if c==141
      title({CFG.fn_data,'Z-over-time'}, 'color','b',...
        'horizontalalignment','left','interpreter','none');
    end
  end
  screen2jpeg([p1,'/',f1,'_z_over_time.jpg'],72)
  pause(2)
  close(gcf)
end

if ~isfield(CFG,'zthres_dvnt')
  zthres_dvnt=5;
else
  zthres_dvnt=CFG.zthres_dvnt;
end
X=permute(cat(3,data.trial{:}),[2 1 3]); % time x chan x trl
IDX1=false(size(X,2),size(X,3));
IDX2=false(size(X,2),size(X,3));
for c=1:size(X,2)
  x=squeeze(X(:,c,:)); % time x trl
  % is the mean value during the trial very different from other trials?
  IDX1(c,:)=abs(zscore(mean(x,1)))>zthres_dvnt;
  % is the value at each timepoint very different from other trials?
  z=zscore(demean(x),[],2);
  IDX2(c,:)=max(abs(z))>zthres_dvnt;
end
isDvnTrl=~~sum(IDX1) | ~~sum(IDX2); % trials including any deviant points
isDvnTrl=isDvnTrl';
if sum(isDvnTrl)
  warning('found %d/%d deviant trials',sum(isDvnTrl),numel(isDvnTrl))
end

if fig
  figure('position',[  -1863          66        1786         881])
  ax=axeslayout([8,20],[0 0 0 0],[0 0 0 0]);
  for c=1:size(X,2)
    axespos(ax,c)
    hold on
    x=squeeze(X(:,c,:)); % time x trl
    z=zscore(denan(x),[],2); % standardized along trials
    plot(data.time{1}, z(:,~isDvnTrl),'k')
    if sum(isDvnTrl)
      plot(data.time{1}, z(:,~~isDvnTrl),'r')
    end
    hold on
    text(0,0,['AG',num2str(c)],'color','w')
    ylim([-10 10])
    if c==141
      title({CFG.fn_data,'Z-over-trial'}, 'color','b',...
        'horizontalalignment','left','interpreter','none');
    end
  end
  screen2jpeg([p1,'/',f1,'_z_over_trial.jpg'],72)
  pause(2)
  close(gcf)
end

end