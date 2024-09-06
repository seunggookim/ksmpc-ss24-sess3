function myfs_viewtkmslices(subid, dname_fssubj, measure, dname_out, caxis, prefix, slice_indices, zoom_factor)
%myfs_viewtkmslices generates surface-overlaid slices
%
% USAGE
%     myfs_viewtkmslices(subid, [dname_fssubj], [measure], [dname_out], [caxis], [prefix], [slice_indices], [zoom_factor])
%
% EXAMPLE
%     % Target directory should be writable and the SUBJECTS_DIR is usually
%     % in /usr/local. So copy the test subject to your home dir:
%     system(['cp -r $SUBJECTS_DIR/bert ~/'])
%
%     % now run
%     myfs_viewtkmslices('bert', getenv('HOME'), 'T1')
%
%     % by default, dname_out = $SUBJECTS/$SUBID/fig/tkmedit/
%
% (cc) 2015, 2017, 2022 dr.seunggoo.kim@gmail.com

if ~exist('dname_out','var')
  dname_out=fullfile(dname_fssubj, subid, 'fig', 'tkmedit');
end
if ~isfolder(dname_out); mkdir(dname_out); end
if ~exist('measure','var')
  measure = 'T1';
end

if ~exist('prefix','var')
  prefix='';
else
  prefix=[prefix,'_'];
end

if ~exist('dname_fssubj','var')
  dname_fssubj=getenv('SUBJECTS_DIR');
else
  setenv('SUBJECTS_DIR',dname_fssubj);
end

if ~exist('slice_indices','var')
  slice_indices=[50:2:200];
else
  if isempty(slice_indices)
    slice_indices=[50:2:200];
  end
end
%
if ~exist('zoom_factor','var')
  zoom_factor=1;
else
  if isempty(zoom_factor)
    zoom_factor=1;
  end
end

slice_indices_padded=[];
for i=1:numel(slice_indices)
  slice_indices_padded = [slice_indices_padded,' ',...
    sprintf('%04i',slice_indices(i))];
end
fname_tcl=fullfile(dname_fssubj,subid,'scripts','tkm_axial.tcl');
fid=fopen(fname_tcl,'w');
fprintf(fid,'# a script for tkmedit \n\n');
fprintf(fid,'LoadVolume %s.mgz\n', measure);
fprintf(fid,'LoadMainSurface 0 %s/%s/surf/lh.pial \n', dname_fssubj, subid);
fprintf(fid,'LoadMainSurface 1 %s/%s/surf/rh.pial \n', dname_fssubj, subid);
fprintf(fid,'foreach SLICE {%s} {\n', slice_indices_padded);
fprintf(fid,'# lh \n');
fprintf(fid,'SetSurfaceLineColor 1 0 0 1 0 \n');
fprintf(fid,'# rh \n');
fprintf(fid,'SetSurfaceLineColor 0 0 0 1 0 \n');
if exist('measure','var')
  switch measure
    case {'uni','T1'}
      fprintf(fid,'SetVolumeMinMax 0 0 400 \n');
    case 'qt1'
      fprintf(fid,'SetVolumeMinMax 0 0 4000 \n');
    case 'bold'
      fprintf(fid,'SetVolumeMinMax 0 0 2000 \n');
    case 'boldresidual'
      fprintf(fid,'SetVolumeMinMax 0 -100 100 \n');
    case 'md'
      fprintf(fid,'SetVolumeMinMax 0 0 0.003 \n');
    case 'uni6k'
      fprintf(fid,'SetVolumeMinMax 0 0 6000 \n');
    case 'b0'
      fprintf(fid,'SetVolumeMinMax 0 0 800 \n');
  end
end
if exist('caxis','var') && ~isempty(caxis)
  fprintf(fid,'SetVolumeMinMax 0 %i %i \n',caxis(1), caxis(2));
end
fprintf(fid,'SetOrientation 1 \n');
fprintf(fid,'SetSlice ${SLICE} \n');
fprintf(fid,'SetZoomLevel %s \n',num2str(zoom_factor));
fprintf(fid,'SetDisplayFlag 22 1 \n');
fprintf(fid,'SetDisplayFlag 3 0 \n');
fprintf(fid,'RedrawScreen \n');
fprintf(fid,'SaveTIFF %s/%s%s_${SLICE}.tif }\n', ...
  dname_out, prefix, subid);
fprintf(fid,'exit \n');
fclose(fid);
system(['tkmedit ',subid,' ',measure,'.mgz -tcl ',fname_tcl]);

end
