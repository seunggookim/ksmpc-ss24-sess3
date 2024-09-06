function fnames_out = myfs_createsurf (subid, subjects_dir, surf2create)
% create additional FreeSurfer surfaces for visualization.
%
% [USAGE]
% fnames_out = myfs_createsurf (subid, [subjects_dir], [cfg])
%
% [INPUT]
% subid        '1xN' subject directory name
% subjects_dir '1xN' directory where you have the subject directory
% surf2create  '1xN' possible inputs are:
%                    'mean'
%                    'semiinflated'
%
% [OUTPUT]


%% CHECK inputs & SET default values
if ~nargin
  help(mfilename)
  return
end

if exist('subjects_dir','var') && ~isempty(subjects_dir)
  setenv('SUBJECTS_DIR',subjects_dir);
else
  subjects_dir = getenv('SUBJECTS_DIR');
end
dname = fullfile(subjects_dir, subid);
if ~isfolder(dname)
  error('%s NOT FOUND!',dname)
end
if ~exist('read_surf','file')
  error('FreeSurfer matlab functions are not in path!')
end


hemis = {'lh.','rh.'};
fnames_out = cellfun(@(x) fullfile(dname,'surf',[x surf2create]), hemis, ...
  'UniformOutput',0);

switch (surf2create)
  case 'mean'
    [surfs] = myfs_readsurfs( subid, subjects_dir, struct('surf',...
      {{'white','pial'}}, 'meas', [], 'annot', []) );
    for ihemi = 1:2
      verts = 0.5*(surfs.white{ihemi}.vertices + surfs.pial{ihemi}.vertices);
      faces = surfs.white{ihemi}.faces;
      write_surf(fnames_out{ihemi}, verts, faces);
    end
    
  case 'semiinflated' % so strange geometry...
    switch (subid)
      case 'fsaverage5'
        niter = 3; dist = 10000;
      case 'fsaverage6'
        niter = 1; dist = 0.1;
      case 'fsaverage'
        niter = 3; dist = 0.1;
      otherwise
        niter = 3; dist = 0.1;
    end
    for ihemi = 1:2
      fname_white = fullfile(dname,'surf',[hemis{ihemi} 'white']);
      system(sprintf('mris_smooth -nw -n %i %s %s',...
        10, fname_white, fnames_out{ihemi}));
      system(sprintf('mris_inflate -no-save-sulc -n %i -dist %i %s %s',...
        niter, dist, fnames_out{ihemi}, fnames_out{ihemi}));
      system(sprintf('mris_smooth -nw -n %i %s %s',...
        10, fnames_out{ihemi}, fnames_out{ihemi}));      
    end
  otherwise
    error('[%s] surf2create=''%s'' NOT RECOGNIZED!',mfilename, surf2create)
end


end