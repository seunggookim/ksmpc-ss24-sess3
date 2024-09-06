function eeg = myeeg_prepmara(fn_eeg, USERUNICA)
% eeg = myeeg_prepmara(fn_eeg, USERUNICA)

fn_out = [fn_eeg(1:end-4),'_ICAMARA.set'];
if isfile(fn_out)
  return
end
eeg = pop_loadset(fn_eeg);


% filtering >1 Hz & downsampling to speed up RUNICA
eeg_hp = pop_eegfiltnew(eeg, 'locutoff',1);
eeg_hp = pop_resample(eeg_hp, 64);
if ~exist('USERUNICA','var') || ~USERUNICA
  % Lazy BINICA creates tempfiles RIGHT IN THE PWD:
  pwd0 = pwd;
  dn = tempname;
  mkdir(dn)
  cd(dn)
  eeg_hp = pop_runica(eeg_hp,'icatype','binica'); % MEMORY efficient!
  cd(pwd0)
else
  eeg_hp = pop_runica(eeg_hp,'icatype','runica'); % but sometimes only this works.
end


% applying the ICA results:
fields = {'icachansind','icawinv','icasphere','icaweights'};
for j = 1:4
  eeg.(fields{j}) = eeg_hp.(fields{j});
end
eeg = eeg_checkset(eeg);


% run mara after BPF [1,30] Hz for artifact detection
eeg = pop_eegfiltnew(eeg,'locutoff',1,'hicutoff',30);
[~,info] = MARA(eeg);
aic = find(info.posterior_artefactprob>0.5);
disp(['MARM found (p>0.5): ',num2str(aic)]);
writematrix(aic',[fn_eeg(1:end-4),'_runica_1-30Hz_MARAaic.txt']);
figure('position',[38    94   884   875],'visible','off')
myeeg_plotcomp(eeg, aic, info.posterior_artefactprob);
export_fig([fn_eeg(1:end-4),'_runica_1-30Hz_MARA.png'],'-r100');
close(gcf);


% subtract the artIC:
eeg = pop_subcomp(eeg, aic);
pop_saveset(eeg, fn_out);

% % Narrower BPF: [1,8] Hz
% eeg = pop_eegfiltnew(eeg,'locutoff',1,'hicutoff',8);
% pop_saveset(eeg, fn_out);
ls(fn_out)
end