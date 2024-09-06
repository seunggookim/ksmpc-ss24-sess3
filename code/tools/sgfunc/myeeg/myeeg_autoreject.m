function EEG_CLEAN = myeeg_autoreject(EEG, keepmeta, suffix)
% EEG_CLEAN = myeeg_autoreject(EEG, [dn_plots], [suffix])

if exist('suffix','var')
  if ~strcmp(suffix(1),'.')
    suffix = ['.',suffix];
  end
end

fn_tmp = [tempname,'.mat'];
save(fn_tmp,'EEG')

[dn_py,~,~] = fileparts(mfilename('fullpath'));
cmd = ['python ',dn_py,'/myar_autorej.py'];
system([cmd,' -i ',fn_tmp])

data = ft_read_data([fn_tmp(1:end-4),'_clean_epo.fif'])*1e6; % volts -> microvolts
meta = load([fn_tmp(1:end-4),'_ar.mat']);

EEG_CLEAN = EEG;
EEG_CLEAN = pop_select(EEG_CLEAN, 'notrial', find(meta.bad_epochs));
assert(EEG_CLEAN.trials == size(data,3))
EEG_CLEAN.data = single(data);
EEG_CLEAN.setname = [EEG_CLEAN.setname, 'cleaned by MNE.autoreject'];
EEG_CLEAN = eeg_checkset(EEG_CLEAN);

if exist('keepmeta','var') && keepmeta
  cmd = ['python ',dn_py,'/myar_cvcurv.py'];
  system([cmd,' -i ',fn_tmp])
  
  dn_trg = [EEG.filepath,'/mne.autoreject',suffix];
  [~,~] = mkdir(dn_trg);
  
  [~,f1,~] = fileparts(EEG.filename);
  prefix_out = [dn_trg,'/',f1];
  fn = dir([fn_tmp(1:end-4),'*png']);
  fn = strcat({fn.folder},'/',{fn.name});
  for j = 1:numel(fn)
    movefile(fn{j}, strrep(fn{j}, fn_tmp(1:end-4), [dn_trg,'/',f1]))
  end
  
  fn = dir([fn_tmp(1:end-4),'_ar*mat']);
  fn = strcat({fn.folder},'/',{fn.name});
  for j = 1:numel(fn)
    movefile(fn{j}, strrep(fn{j}, fn_tmp(1:end-4), [dn_trg,'/',f1]))
  end
  
  %% more sanity check figures
  figure('visible','off')
  imagesc(meta.labels_epochs_by_channels')
  hold on
  scatter(find(meta.bad_epochs), ...
    EEG.nbchan*ones(sum(meta.bad_epochs))+0.5,'ro')
  axis ij
  xlabel('Trials'); ylabel('Chan');
  colormap([1 1 1; 0 1 0; 1 0 0])
  exportgraphics(gcf,[dn_trg,'/',f1,'_labels.png'],'Resolution',150)
  close(gcf)

  n = ceil(sqrt(EEG.trials));
  ax = axeslayout([1 1]*n);
  figure('position',[1 1 100*n 100*n],'visible','off')
  for itrl = 1:EEG.trials
    axespos(ax, itrl)
    hold on
    plot(EEG.times, EEG.data(:,:,itrl)','k');
    
    remaining_trials = find(~meta.bad_epochs);
    if ismember(itrl, remaining_trials)
      jtrl = find(remaining_trials==itrl);
      plot(EEG_CLEAN.times, EEG_CLEAN.data(:,:,jtrl)','color','r');
    else
      text(0,0,'rejected','color','r')
    end
    title(sprintf('Trial #%02i',itrl))
    ylim([-100 100])
    if itrl == 1
      ylabel('\muV'); xlabel('ms')
    end
  end
  exportgraphics(gcf,[dn_trg,'/',f1,'_trials.png'],'Resolution',150)
  close(gcf)
  
end

end
