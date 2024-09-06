function myft_plot_timelockrms_eeg(timelock, titlestr)
%%
clf;
set(gcf,'position',[  -926   532   635   398])
set(0,'DefaultAxesColorOrder',brewermap(9,'set1'));
innermargin=[.1 .05 .1 .13];
ax1=axeslayout([1 1],innermargin,[0 0 0 0]);
axespos(ax1,1)
plot(timelock.time, timelock.avg,'color',[0 0 0 0.5])
ylabel('Potential [uV]')
xlabel('Time [s]');
ylim0=[min(timelock.avg(:)) max(timelock.avg(:))];
ylim([ylim0(1) ylim0(2)*1.8])
grid on
c0=max(abs(timelock.avg(:)));
c0=[-c0 c0]*1;
if exist('titlestr','var')
  title(titlestr,'interpreter','none','fontsize',9)
end

% compute rms
y=rms(timelock.avg);
dt=timelock.time(2)-timelock.time(1);
[pks,idx_pks]=findpeaks(y,'MinPeakDistance',0.050/dt);
idx_idx=timelock.time(idx_pks)<0;
idx_pks(idx_idx)=[];
pks(idx_idx)=[];
t_pks=timelock.time(idx_pks);
hold on;
plot(timelock.time,y,'color',[.9 0 0],'linewidth',0.5)

% connect peaks & topoplots
nPks=numel(idx_pks);
y0=ylim;
for j=1:nPks
  hold on
  topo_x0=timelock.time(1)+range(timelock.time)/nPks*(j-1+0.5);
  topo_y0=y0(2)*0.75;
  %line([t_pks(j) topo_x0]',[pks(j) topo_y0]', 'color',[.2 .2 0],'linewidth',0.5)
  line([t_pks(j) topo_x0]',[pks(j) topo_y0]', 'color',[.2 0 .2],'linewidth',0.5)
end

ax2=axeslayout([3 nPks],[0 0 0.1 0.15],innermargin);
cfg=struct('highlightcolor','k','zlim',c0);
cfg.layout=ft_prepare_layout(struct('layout','easycapM24.mat'));
for j=1:nPks
  axespos(ax2,j)
  cfg.xlim=[t_pks(j)-dt t_pks(j)+dt];
  myft_topo(cfg, timelock)
end

end