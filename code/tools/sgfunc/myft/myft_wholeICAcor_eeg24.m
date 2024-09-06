function myft_wholeICAcor_eeg24(fn_data,CFG)
[p1,f1,e1] = myfileparts(fn_data);
if ~exist('CFG','var')
  CFG=[];
end
if ~isfield(CFG,'bpfreq')
  CFG.bpfreq = [0.1 50]; % online filtering was bandpass=[1-200] notch=[60]
end
freqsuffix = [num2str(CFG.bpfreq(1)),'-',num2str(CFG.bpfreq(2)),'Hz'];
prefix = [p1,'/',f1,'_',freqsuffix];
fn_ds = [prefix,'_ds.mat'];
if ~exist(fn_ds,'file')
  %% read data with a given filtering
  cfg=[];
  cfg.dataset = fn_data;
  cfg.channel = 1:24;
  cfg.demean = 'yes'; % demean by baseline
  cfg.baselinewindow = [-Inf Inf]; % from beginning to end
  if CFG.bpfreq(1)==0 || CFG.bpfreq(2)==inf
    error('now write down high-pass or low-pass only config HERE!')
  end
  cfg.bpfilter = 'yes';
  cfg.bpfreq=CFG.bpfreq;
  if CFG.bpfreq(1)<1
    cfg.bpfilttype='firws';  
  else
    cfg.bpfilttype='but';
    cfg.bpfiltord=4;
  end
  cfg.bpfiltdir='twopass';
  cfg.precision='single';
  data = ft_preprocessing(cfg);
  % downsampling
  cfg = [];
  cfg.resamplefs = 200;
  data = ft_resampledata(cfg, data);
  save(fn_ds,'data')
else
  load(fn_ds,'data')
end
%% ICA
fn_comp=[prefix,'_ICs.mat'];
nCmp = 23;
if ~exist(fn_comp,'file')
  cfg=[];
  cfg.method = 'runica';
  cfg.numcomponent = nCmp;
  cfg.channel = 1:24;
  comp = ft_componentanalysis(cfg, data);
  save(fn_comp,'comp')
else
  load(fn_comp,'comp')
end
%% find candidates
fn_png=[prefix,'_ICs.png'];
if ~exist(fn_png,'file')
  figure(1);clf;
  set(1,'position',[ -1909           1        1890         379]);
  ax=axeslayout([2 20],[0 0 0 0],[0 0 0 0]);
  ft_info off
  lay=ft_prepare_layout(struct('layout','easycapM24.mat'));
  for c=1:nCmp
    axespos(ax,c)
    cfg=[];
    cfg.component=c;
    cfg.layout=lay;
    cfg.comment='no';
    cfg.markersymbol='.';
    cfg.title='off';
    ft_topoplotIC(cfg, comp)
    text(0,0.8,num2str(c),...
      'horizontalAlignment','center','fontweight','bold');
  end
  ft_info on
  axes('position',[0 0 1 .5])
  title(fn_comp,'Interpreter','none')
  axis off
  colormap(flipud(brewermap(256,'spectral')))
  screen2png(fn_png,150)
else
  figure(1);clf;
  axes('position',[0 0 1 1])
  image(imread(fn_png));
  axis image; axis off;
  set(1,'position',[ -1909           1        1890         379]);
end
if isfield(CFG,'runICAonly')&&CFG.runICAonly
  return
end
%% Now interactive session:
alldone=0;
while ~alldone
  comptype={'EOG','ECG'};
  Idx_artcomp={};
  for k=1:2
    %fn_tsv=[p1,'/',f1,'_',comptype{k},'comp.tsv'];
    fn_tsv=[prefix,'_',comptype{k},'comp.tsv'];
    go=0;
    while ~go
      idx_artcomp = input(['[?] ',comptype{k},'-components (e.g. [2 3])? ']);
      if prod(~isnan(idx_artcomp))
        HF=[];
        for j=1:numel(idx_artcomp)
          HF=[HF figure('position',[500*(j-1)-1900 650 400 300])];
          myft_checkIC1(comp,idx_artcomp(j));
        end
        go=input('[?] Good to go (1/0)? ');
        close(HF)
      else
        go=1;
      end
    end
    dlmwrite(fn_tsv,idx_artcomp,'delimiter','\t');
    Idx_artcomp{k}=idx_artcomp;
  end
  
  %% detect a spike from each timeseries:
  % for EOG:
  cfg_preset{1}=struct('hilbert','yes', 'bpfilter','yes', ...
    'bpfilttype','but','bpfreq',[1 10], 'bpfilord',4, 'zcutoff',5);
  % for ECG:
  cfg_preset{2}=struct('hilbert','yes', 'hpfilter','yes', ...
    'hpfilttype','firws','hpfreq',[1], 'zcutoff',3);
  Events={};
  for k=1:2 % k=1 EOG, k=2 ECG
    %unix(['rm ',p1,'/',f1,'_',comptype{k},'comp*png']); % just reset
    unix(['rm ',prefix,'_',comptype{k},'comp*png']); % just reset
    idx_artcomp=Idx_artcomp{k};
    if prod(~isnan(idx_artcomp))
      events=[];
      for j=1:numel(idx_artcomp)
        c=idx_artcomp(j);
        cfg=cfg_preset{k};
        cfg.channel=c;
        dat = ft_preprocessing(cfg, comp);
        %%
        figure
        set(gcf,'position',[ -1843         456        1287         409])
        ax2=axeslayout([2 1],[.02 .02 .15 .1]);
        axespos(ax2,1);
        [pks,idx]=findpeaks(zscore(dat.trial{1}),...
          'MinPeakHeight',cfg.zcutoff,...
          'MinPeakDistance',dat.fsample*0.5);
        t=dat.time{1}(idx);
        hold on;
        plot(dat.time{1}, zscore(dat.trial{1}),'k');
        scatter(t,pks,'v','markerFaceColor','r')
        xlabel('Time [s]'); ylabel('|Z|');
        title({fn_data,[comptype{k},'-suspected: ',...
          'found ',num2str(numel(t)),' peaks, ',...
          'dist=',num2str(mean(diff(t))),' \pm ',...
          num2str(std(diff(t))),' sec']}, 'Interpreter','none');
        grid on
        xlim([dat.time{1}(1) dat.time{1}(end)])
        
        ax=axeslayout([2 5],[.1 .1 .2 .2]);
        axespos(ax,6);
        ft_topoplotIC(struct('component',c,'layout','easycapM24.mat',...
          'comment','no','markersymbol','.'), comp);
        colormap(gca,flipud(brewermap(256,'spectral')))
        
        % re-epoch IC timeseries
        trl=round(comp.fsample.*[t'-1 t'+1 0*t'-1]);
        trl(trl(:,1)<-trl(1,3),:)=[];
        trl(trl(:,2)>(-trl(1,3)+size(comp.time{1},2)),:)=[];
        if size(trl,1)>200
          warning(['Too many trials: randomly select 200']);
          idx=randperm(size(trl,1));
          trl=trl(idx(1:200),:);
        end
        cfg=[];
        cfg.trl=trl;
        epoch = ft_redefinetrial(cfg, comp);
        
        axespos(ax,7);
        [PSD,f]=pwelch(comp.trial{1}(c,:),1024,[],[],comp.fsample);
        h=plot(f, PSD, 'b');
        xlabel('Frequency [Hz]'); ylabel('W-PSD [au^2/Hz]')
        set(gca,'xscale','log','yscale','lin','XMinorGrid','on',...
          'Yminorgrid','off')
        
        axespos(ax,8)
        E = cat(3,epoch.trial{:});
        errorplot(epoch.time{1}, denan(squeeze(E(c,:,:)))');
        xlabel(['Peak-related time [s]']); ylabel('IC weight [au]')
        grid on;
        
        axespos(ax,9)
        mypcolor(epoch.time{1},1:size(trl,1),denan(squeeze(E(c,:,:)))', ...
          struct('grid',0))
        xlabel(['Peak-related time [s]']); title('IC weight [au]')
        ylabel('Trial index'); set(gca,'ydir','rev')
        c0=caxis;
        caxis([-max(abs(c0)) max(abs(c0))]*0.6)
        colormap(gca,bipolar)
        
        axespos(ax,10)
        hold on;
        cfg=[];
        cfg.trl=trl;
        epoch = ft_redefinetrial(cfg, data);
        E = cat(3,epoch.trial{:});
        plot(epoch.time{1}, squeeze(nanmean(E(1:157,:,:),3)),'color','k');
        
        cfg=[];
        cfg.component=c;
        data_recon = ft_rejectcomponent(cfg, comp, data);
        cfg=[];
        cfg.trl=trl;
        epoch = ft_redefinetrial(cfg, data_recon);
        E=cat(3,epoch.trial{:});
        plot(epoch.time{1}, squeeze(nanmean(E(1:157,:,:),3)),'color','r');
        xlabel(['Peak-related time [s]']); ylabel('Potential [uV]')
        grid on;
        title(['Before/{\color[rgb]{1 .2 .2}after} IC rejection'])
        fn_png=[prefix,...
          '_',comptype{k},'comp',num2str(idx_artcomp(j)),'.png'];
        screen2png(fn_png,150)
        close(gcf)
        events=[events t];
      end
    else
      events=nan;
    end
    Events{k}=unique(events);
    fn_tsv=[prefix,'_',comptype{k},'_events.tsv'];
    dlmwrite(fn_tsv,Events{k},'delimiter','\t');
    ls(fn_tsv)
  end
  alldone=input('[?] All done(0/1)? ');
end
close(1)
%% Reject component:
cfg=[];
cfg.component = [Idx_artcomp{1}, Idx_artcomp{2}];
data = ft_rejectcomponent(cfg, comp, data);
data.grad.balance.current='none';
fn_datacor=[prefix,'_ICAcor.mat'];
save(fn_datacor,'data')
ls(fn_datacor)

end

%% SUBFUNCTIONS: 
function myft_checkIC1(comp,c)
%%
nTrl=numel(comp.trial);

ax=axeslayout([2 3]);
axespos(ax,1)
cfg=[];
cfg.component=c;
cfg.layout='easycapM24.mat';
cfg.comment='no';
cfg.markersymbol='.';
ft_topoplotIC(cfg, comp)
colormap(gca,flipud(brewermap(256,'spectral')))

if nTrl>1
  
  X=cat(3,comp.trial{:});
  axespos(ax,2)
  hold on
  mypcolor(comp.time{1},1:numel(comp.trial),squeeze(X(c,:,:))');
  caxis0=caxis;
  caxis([-max(abs(caxis0)) max(abs(caxis0))]./2);
  %colormap(gca,jet)
  colormap(gca,flipud(brewermap(256,'RdYlBu')))
  %colormap(gca,sgcolormap('BCWYR256'));
  ylabel('Trial index');
  xlabel('Time [s]');
  title('IC weights [AU]')
  
  axespos(ax,5)
  hold on
  errorplot(comp.time{1}, squeeze(X(c,:,:))')
  xlabel('Time [s]'); ylabel('IC wegith [AU]');
  grid on; box on;
  
  axespos(ax,3)
  scatter(squeeze(mean(X(c,:,:),2)),1:nTrl,'filled',...
    'markerfacecolor',[0 0 0],'markerFaceAlpha',0.2)
  xlabel('Mean weight [AU]'); ylabel('Trial index')
  grid on; box on;
  
  axespos(ax,6)
  scatter(squeeze(std(X(c,:,:),[],2)),1:nTrl,'filled',...
    'markerfacecolor',[0 0 0],'markerFaceAlpha',0.2)
  xlabel('Std weight [AU]'); ylabel('Trial index')
  grid on; box on;
  
  axespos(ax,4)
  [~,f]=pwelch(comp.time{1}(1,:),[],[],[],comp.fsample);
  PSD=zeros(nTrl,numel(f));
  for t=1:nTrl % 141 samples
    [PSD(t,:),f]=pwelch(comp.trial{t}(c,:),[],[],[],comp.fsample);
  end
  [pks,idx]=findpeaks(mean(PSD));
  [~,idx_idx]=max(pks);
  idx=idx(idx_idx);
  pk=pks(idx_idx);
  h=errorplot(f, PSD);
  xlabel('Frequency [Hz]'); ylabel('W-PSD [au^2/Hz]')
  set(gca,'xscale','log','yscale','lin','XMinorGrid','on','Yminorgrid','off')
  hold on;
  scatter(f(idx),pk,'bx')
else
  axespos(ax,2)
  [PSD(1,:),f]=pwelch(comp.trial{1}(c,:),[],[],[],comp.fsample);
  [pks,idx]=findpeaks(PSD);
  [~,idx_idx]=max(pks);
  idx=idx(idx_idx);
  pk=pks(idx_idx);
  h=errorplot(f, PSD);
  xlabel('Frequency [Hz]'); ylabel('W-PSD [au^2/Hz]')
  set(gca,'xscale','log','yscale','lin','XMinorGrid','on','Yminorgrid','off')
  hold on;
  scatter(f(idx),pk,'bx')
  
  ax2=axeslayout([2 1],[.1 .05 .1 .15],[.02, .3, .02, .02]);
  axespos(ax2,2)
  x=comp.trial{1}(c,:);
  plot(comp.time{1}, x)
  xlabel('Time [s]');ylabel('IC weights [AU]');
  title('Whole data')
  
  axespos(ax,6)
  [~,s]=max(x(401:end-400));
  s=s+400;
  plot(comp.time{1}, x)
  xlim([comp.time{1}(s-400) comp.time{1}(s+400)])
  xlabel('Time [s]');ylabel('IC weights [AU]');
  title('Zoomed around max')
  
  axespos(ax,3)
  [acor,lag]=xcorr(x,4*comp.fsample);
  plot(lag/comp.fsample,acor,'b')
  xlim([-4 4])
  xlabel('Lag [s]'); ylabel('Autocorr');
end

end