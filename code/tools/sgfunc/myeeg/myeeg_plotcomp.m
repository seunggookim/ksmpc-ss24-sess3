function myeeg_plotcomp(EEG, aic, prob)
if ~exist('aic','var')
  aic = [];
end

ncomp = size(EEG.icawinv,2);
k = ceil((ncomp+3)^0.5);
ax = axeslayout([k k],'tight','tight');
for icomp = 1:ncomp
  axespos(ax,icomp)
  topoplot(EEG.icawinv(:, icomp), EEG.chanlocs);
  h = title(['IC',num2str(icomp)]);
  if ismember(icomp, aic)
    h.BackgroundColor = 'r';
    h.Color = 'w';
  end
  caxis([-5 5])
end
set(gcf,'color','w')
colormap(flipud(brewermap(256,'spectral')))

axespos(ax,ncomp+1)
XtX = EEG.data(1:EEG.nbchan,:)*EEG.data(1:EEG.nbchan,:)';
[~,D] = eig(XtX);
plot(cumsum(diag(D))/sum(diag(D)),'bo-')
grid on; ylim0 = ylim;
text(ncomp*0.95, ylim0(1), ...
  {'X''X',['rank=',num2str(rank(XtX))]},...
  'horizontalAlignment','right', 'verticalAlignment','bottom',...
  'fontsize',7)

axespos(ax,ncomp+2)
if isempty(EEG.icaact)
%   EEG.icaact = EEG.icaweights * EEG.icasphere * EEG.data(1:size(EEG.icaweights,1),:,:);
  EEG.icaact = eeg_getdatact(EEG);
end
d = size(EEG.icaact);
if numel(d)<3
  d(3) = 1 ;
end
A = reshape(EEG.icaact,[d(1), d(2)*d(3)]);
AtA = A*A';
[~,D] = eig(AtA);
plot(cumsum(diag(D))/sum(diag(D)),'ro-')
grid on; ylim0 = ylim;
text(ncomp*0.95, ylim0(1), ...
  {'A''A',['rank=',num2str(rank(AtA))]},...
  'horizontalAlignment','right', 'verticalAlignment','bottom',...
  'fontsize',7)

if exist('prob','var')
  axespos(ax,ncomp+3)
  plot(sort(prob),'go-')
  grid on; ylim0 = ylim;
  text(ncomp*0.95, ylim0(1), ...
  'Pr(art|feat)',...
  'horizontalAlignment','right', 'verticalAlignment','bottom',...
  'fontsize',7)
end

end