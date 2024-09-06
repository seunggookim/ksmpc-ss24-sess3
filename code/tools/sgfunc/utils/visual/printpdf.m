function printpdf(fname_pdf, dpi)
% printpdf(fname_pdf, [dpi])
%

drawnow;
fig = gcf;
fig.PaperPositionMode = 'auto';
fig.InvertHardcopy = 'off';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
fig.Renderer = 'Opengl';
if ~exist('dpi','var')
  dpi = 600;
end
print(fig,'-dpdf',['-r',num2str(dpi)],'-bestfit', fname_pdf)
try 
  unix(['pdfcrop ',fname_pdf,' ',fname_pdf])
catch
  warning('pdfcrop cannot be run')
end
end