function gspdf2png(FnPdf, Dpi)
% gspdf2png(FnPdf, Dpi)

if not(isunix)
  error('Use UNIX!')
end
if not(exist('Dpi','var')), Dpi=300; end

system(sprintf('gs -sDEVICE=pngalpha -r%i -o %s %s', Dpi, strrep(FnPdf,'.pdf','.png'), FnPdf));
end