function myeeg_topocomp(EEG)
% myeeg_topocomp(EEG)
% (cc) sgKIM, 2019.
figure('position',[1921          78        1470         896])
ax   = axeslayout([6 11],'tight','tight');
for icomp=1:63
  axespos(ax,icomp)
  topoplot(EEG.icawinv(:, icomp), EEG.chanlocs);
  title(['IC',num2str(icomp)]);
end
set(gcf,'color','w')
colormap(flipud(brewermap(256,'spectral')))

end