function myft_plot_timelock(timelock, titlestr)
%%
clf;
set(gcf,'position',[  -926   532   635   398])
set(0,'DefaultAxesColorOrder',brewermap(9,'set1'));
innermargin=[.1 .05 .1 .13];
ax1=axeslayout([1 1],innermargin,[0 0 0 0]);
axespos(ax1,1)
plot(timelock.time, timelock.avg)
ylabel(['Field [',timelock.grad.chanunit{1},']'])
xlabel('Time [s]');
ylim0=ylim;
ylim([ylim0(1) ylim0(2)*1.5])
grid on
c0=max(abs(timelock.avg(:)));
c0=[-c0 c0]*0.2;
if exist('titlestr','var')
  title(titlestr,'interpreter','none','fontsize',9)
end

%layout=ft_prepare_layout(struct('layout','yokogawa160_helmet.mat'));
layout=ft_prepare_layout(struct('grad',timelock.grad));
T=get(gca,'xtick');
ax2=axeslayout([1 numel(T)-1],[0 0 0 0.7],innermargin);
for i=1:(numel(T)-1)
  axespos(ax2,i)
  cfg=struct('layout',layout, 'comment','no', 'title','off', ...
    'markersymbol','.', 'zlim',c0);
  cfg.xlim=[T(i) T(i+1)];
  ft_warning off; ft_notice off; ft_info off
  ft_topoplotER(cfg, timelock);
  %axis([-.9 .9 -.8251 .71])
end


end