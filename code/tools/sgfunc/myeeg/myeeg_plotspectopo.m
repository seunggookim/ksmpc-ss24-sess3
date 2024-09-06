function fn_png = myeeg_plotspectopo(eeg)
figure('visible','off')
pop_spectopo(eeg, 1, [0  eeg.pnts], 'EEG' , ...
  'percent', 50, 'freq', [6 12 25], 'freqrange',[2 30],'electrodes','on');
colormap(flipud(brewermap(256,'RdBu')));
h = get(gcf,'children');
set(h(2),'color','w')
hh = get(h(1),'children');
hh.CDataMapping = 'scaled';
set(gcf,'color','w')
h(6).FontSize = 10;
[~,f1,~] = fileparts(eeg.filename);
fn_png = fullfile(eeg.filepath, [f1,'_spectopo.png']);
export_fig(fn_png,'-r120')
close(gcf)
end
