function myeeg_ploticlabel(eeg, cfg)

ICLabel = eeg.etc.ic_classification.ICLabel;
fn = {};

ncomp = size(eeg.icaweights,1);
assert(ncomp == size(ICLabel.classifications,1),'#comp not match!')

% Topoplot + IC label probabilities
figure('position',[38    63   950   881],'visible','off')
ax = axeslayout([1 1]*ceil(sqrt(ncomp*2)),'tight',[.02, .02, .03, .01]);
for icomp = 1:ncomp
  axespos(ax,2*(icomp-1)+1)
  topoplot(eeg.icawinv(:, icomp), eeg.chanlocs, 'electrodes','off');
  h = title(['IC',num2str(icomp)]);
  caxis([-5 5])
  
  axespos(ax,2*(icomp))
  h = bar(ICLabel.classifications(icomp,:));
  h.CData = brewermap(7,'spectral');
  h.FaceColor = 'flat';
  ylim([0 1]);
  set(gca,'xtick',1:7, 'xticklabel',cellfun(@(x) x(1), ICLabel.classes, ...
    'uni',0))
  axis square
end
set(gcf,'color','w')
axes('position',[0 0 1 0.98]); axis off;
if isfield(cfg,'title')
  title(cfg.title,'interp','none','fontsize',10)
else
  title(fullfile(eeg.filepath,eeg.filename),'interp','none','fontsize',10)
end
colormap(flipud(brewermap(256,'spectral')))

fn{1} = [tempname,'.png'];
exportgraphics(gcf, fn{1},'Resolution',150)
close(gcf)

% ERP images
epoch = pop_resample(eeg,250);
epoch = pop_eegfiltnew(epoch,[0.5 20]);
events = unique({eeg.event(contains({eeg.event.code},'Stimulus')).type});
epoch = pop_epoch(epoch,events,cfg.toi_s);
erp = pop_rmbase(epoch,[cfg.toi_s(1) 0]);
erp = myeeg_icaact(erp);
erp.times = erp.times./1000;
%%
figure('position',[938   352   316   521],'visible','off')
ax = axeslayout([8 2],[.15, .0, .2, .2],[.02, .02, .0, .02]);
axespos(ax,1)
ichan = eeg_chaninds(eeg, 'FCz');
imagesc(erp.times, 1:erp.trials, squeeze(erp.data(ichan,:,:))')
ylabel('Trial')
colorbar
caxis([-100 100])
colormap(gca, flipud(brewermap(128,'RdBu')))
title('FCz','fontsize',6)
set(gca,'fontsize',6)

axespos(ax,2)
ichans = find(contains({eeg.chanlocs.labels}, 'EOG'));
eognorm = squeeze(vecnorm(erp.data(ichans,:,:),2,1));
imagesc(erp.times, 1:erp.trials, eognorm')
ylabel('Trial')
colorbar
caxis([0 300])
colormap(gca, flipud(gray))
title('EOG-l2norm','fontsize',6)
set(gca,'fontsize',6)

[~,idxs] = sort(sum(ICLabel.classifications(:,2:6)), 'descend');
idxs = idxs + 1;
for i = 1:2
  [~,idx] = sort(ICLabel.classifications(:,idxs(i)),'descend');
  for j = 1:6
%     subplot(8,2,2+2*(j-1)+i)
    axespos(ax,2+2*(j-1)+i)
    imagesc(erp.times, 1:erp.trials, squeeze(erp.icaact(idx(j),:,:))')
    x = ICLabel.classes{idxs(i)};
    title(sprintf('[IC%i] %s=%.2f',idx(j), ...
      x(1:min([4 numel(x)])),...
      ICLabel.classifications(idx(j),idxs(i))), 'fontsize',6)
    colormap(gca,flipud(brewermap(128,'spectral')))
    colorbar
    set(gca,'fontsize',6)
  end
  xlabel('Time [s]')
end

ax = axeslayout([10 1],'tight');
axespos(ax,10)
h = bar(sum(eeg.etc.ic_classification.ICLabel.classifications(:,1:7)>0.6));
h.CData = brewermap(7,'spectral');
h.FaceColor = 'flat';
set(gca,'xtick',1:7, 'xticklabel',cellfun(@(x) x(1), ICLabel.classes, ...
  'uni',0))
ylabel('# comp > 60%')
set(gca,'fontsize',6)

set(gcf,'color','w')

%%
fn{2} = [tempname,'.png'];
exportgraphics(gcf, fn{2},'resolution','250')
close(gcf)

[~,f1,~] = fileparts(eeg.filename);
fn_out = fullfile(eeg.filepath, [f1,'_iclabel.png']);
imageconcat(fn, fn_out, 2)
end