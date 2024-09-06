function [fname_out] = savevtc2nifti(data, InfoVTC, fname_ref, fname_out)
% [fname_out] = savevtc2nifti(data, InfoVTC, fname_ref, fname_out)
compressed = contains(fname_out,'.gz');

map = InfoVTC.voxVTC*0;
map(InfoVTC.voxVTC) = data;
map = reshape(map, InfoVTC.DimVTC);

if isstruct(fname_ref)
  info = fname_ref;
else
  info = niftiinfo(fname_ref);
end
info.Filename = fname_out;
info.Datatype = class(map);
info.Description = 'VTCdata';

[p1,f1,~] = fileparts(strrep(fname_out,'.gz',''));
niftiwrite(map, fullfile(p1,f1), info, 'Compressed', compressed)

end
