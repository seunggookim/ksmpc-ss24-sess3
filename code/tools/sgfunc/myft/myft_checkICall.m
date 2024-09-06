function myft_checkICall(comp)
%{
Basic cues are:

1. ECG topographies (usually in pairs) are very smooth
2. ECG peaks are somewhat regular
3. EOG topography peaks in the occular channels
4. Eye blink should be very sharp and strong

%}
dp=get(0,'DefaultFigurePosition');
figure('position',[dp(1)   43        1380         912]);
% % A. find spikes
% X=cat(2,comp.trial{:});
% subplot(5,5,1)
nTrl=numel(comp.trial);
nCmp=size(comp.topo,2);
% t_alltrl=[];
% for t=1:nTrl
%   t_alltrl=[t_alltrl comp.time{1}+(t-1)*comp.time{1}(end)];
% end
% plot(t_alltrl,abs(zscore(X')))
% xlabel('Time [s]: all epoches concat')
% ylabel('|Z|')
% idx_cmp=find(max(abs(zscore(X')))>10);
% title(['|Z|>10: ',num2str(idx_cmp)])
% for k=1:numel(idx_cmp)
%   subplot(5,5,1+k)
%   cfg=[];
%   cfg.component=idx_cmp(k);
%   cfg.grad=comp.grad;
%   cfg.comment='no';
%   cfg.markersymbol='.';
%   ft_topoplotIC(cfg, comp)
%   colormap(gca,flipud(brewermap(256,'spectral')))
% end

% B. find peaks in spectrogram
subplot(5,5,1);
hold on;
[psd,f]=pwelch(comp.time{1}(1,:),[],[],[],comp.fsample);
psd0=psd*0;
peakfreq=zeros(1,nCmp);
for c=1:nCmp
  psd=psd0;
  for t=1:nTrl % 141 samples
    [psd1,f]=pwelch(comp.trial{t}(c,:),[],[],[],comp.fsample);
    psd=psd+psd1;
  end
  psd=psd./nTrl;
  [pks,idx]=findpeaks(psd);
  [~,idx_idx]=max(pks);
  idx=idx(idx_idx);
  peakfreq(c)=f(idx);
  plot(f, psd);
end
xlabel('Frequency [Hz]'); ylabel('W-PSD [au^2/Hz]')
set(gca,'xscale','log','yscale','lin','XMinorGrid','on','Yminorgrid','off')

subplot(5,5,2)
hist(peakfreq,nCmp)
xlabel('Peak frequency [Hz]')
ylabel('# of compo')
idx_cmp=find(peakfreq<=1); % check out these slow compos
title(['Peak at <=1 Hz: ',num2str(idx_cmp)])
for k=1:numel(idx_cmp)
  subplot(5,5,2+k)
  cfg=[];
  cfg.component=idx_cmp(k);
  cfg.grad=comp.grad;
  cfg.comment='no';
  cfg.markersymbol='.';
  ft_topoplotIC(cfg, comp)
  colormap(gca,flipud(brewermap(256,'spectral')))
end

end