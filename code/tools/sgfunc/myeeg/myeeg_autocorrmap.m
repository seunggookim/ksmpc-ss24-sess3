function EEGout = myeeg_autocorrmap (fname_in, cfg)
% EEGout = myeeg_autocorrmap (fname_in, cfg)
% 
% fname_in   - EEG file
% post_thres - posterior threshold (default = 0.5)

EEG = pop_loadset(fname_in);
[p1,f1,e1] = fileparts(fname_in);
fname_out = fullfile(p1,[f1,'_cor',num2str(post_thres),e1]);
fname_mat = fullfile(p1,[f1,'_cor.mat']);
% givewarning
if exist(fname_out,'file')
  ls(fname_out)
  if nargout
    EEGout = pop_loadset(fname_out);
  end
  return
end

%% Run CORRMAP
cfg0 = [];
% defintion of template datasets and components for each of the three
% components that should be used by the semi-automatic corrmap algorithm
cfg0.blink_set = 8;    % identified set for the template eye blink component
cfg0.blink_comp = 1;   % identified template component included in eye_set
cfg0.heart_set = 8;    % identified set for the template heart beat  component
cfg0.heart_comp = 10;  % identified template component included in heart_set
cfg0.eyemov_set = 8;   % identified set for the template eye movement component
cfg0.eyemov_comp = 17; % identified template component included in eye_set
% Parameters are from [https://doi.org/10.3389/fnins.2018.00309]




if ~exist(fname_mat,'file')
  [~,info] = MARA(EEG);
  save(fname_mat, 'info')
else
  load(fname_mat, 'info')
end
post = info.posterior_artefactprob;

%% Creating component plots
fname_png = fullfile(p1,[f1,'_MARA',num2str(post_thres),'.png']);
if ~exist(fname_png, 'file')
  figure('position',[1921, 78, 1470, 896], 'visible','off')
  ax = axeslayout([6 11],'tight','tight');
  ncomp = size(EEG.icawinv,2);
  nchan = EEG.nbchan;
  for icomp = 1:ncomp
    axespos(ax,icomp)
    topoplot(EEG.icawinv(:, icomp), EEG.chanlocs);
    h = title(sprintf('IC%02i [%.2f]', icomp, post(icomp)));
    if post(icomp) > post_thres
      h.BackgroundColor = 'r';
      h.Color = 'w';
    end
    caxis([-5 5])
  end
  set(gcf,'color','w')
  colormap(flipud(brewermap(256,'spectral')))
  
  axespos(ax,ncomp+1)
  [~,D] = eig(EEG.data(1:nchan,:)*EEG.data(1:nchan,:)');
  eigval = sort(diag(D),'descend');
  plot(cumsum(eigval)/sum(eigval),'bo-')
  grid on; ylim0 = ylim;
  text(ncomp*0.95, ylim0(1), {'Cumulative','normalized','eigenvalues of X''X',...
    'over eigvec',['rank=',num2str(rank(EEG.data(1:nchan,:)*EEG.data(1:nchan,:)'))]},...
    'fontsize',10,'horizontalAlignment','right', 'verticalAlignment','bottom')
  
  axespos(ax,ncomp+2)
  if isempty(EEG.icaact)
    EEG.icaact = EEG.icaweights * EEG.icasphere * EEG.data(1:ncomp,:,:);
  end
  [~,D] = eig(EEG.icaact(1:ncomp,:)*EEG.icaact(1:ncomp,:)');
  eigval = sort(diag(D),'descend');
  plot(cumsum(eigval)/sum(eigval),'ro-')
  grid on; ylim0 = ylim;
  text(ncomp*0.95, ylim0(1), {'Cumulative','normalized','eigenvalues of A''A',...
    'over eigvec',['rank=',num2str(rank(EEG.icaact(1:ncomp,:)*EEG.data(1:ncomp,:)'))]},...
    'fontsize',10,'horizontalAlignment','right', 'verticalAlignment','bottom')
  
  axespos(ax,ncomp+3)
  hold on;
  sortedpost = sort(post);
  plot(sortedpost, 'b.','linestyle','none')
  idx = find(sortedpost > post_thres);
  plot(idx, sortedpost(idx), 'r.', 'linestyle','none')
  text(2,1,{'Sorted','posterior','over comp'}, 'verticalAlignment','top',...
    'fontsize',10,'interpret','none')
  line([0 nchan+1]', post_thres*[1 1]', 'color','r')
  line(sum(post < post_thres)*[1 1]', [0 1]', ...
    'color',[.3 .3 .3],'linestyle','-','linewidth',0.5)
  axis([1 nchan 0 1])
  box on
  
  export_fig(fname_png)
end

%%
artcomps = find(post > post_thres);
disp(['Removing component(s) with posterior > ',num2str(post_thres),': ',num2str(artcomps)])
EEGout = pop_subcomp(EEG, artcomps);
pop_saveset(EEGout, fname_out);
fprintf('Saving in ')
ls(fname_out)
if ~nargout
  clear EEGout
end
% givewarning
end

% function givewarning
% warning('***THIS IS ONLY FOR PRELIMINARY ANALYSIS!!!!***')
% warning('***YOU HAVE TO INSPECT ALL COMPONENTS CAREFULLY FOR ACTUAL ANALYSIS!!!***')
% end
