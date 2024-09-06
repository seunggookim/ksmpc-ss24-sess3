function [H, cfg] = myfs_view (spaces, data, cfg)
%MYFS_VIEW visualies scalar data on surfaces (and volume)
%
% [USAGE]
% myfs_view (spaces, data)
% H = myfs_view (spaces, data)
% [H, cfg] = myfs_view (spaces, data, cfg)
%
%
% [H, cfg] = myfs_view ([], NiftiFileNameInMNISpace)
% automatically view assuming the NIFTI volume is in the MNI space.
%
% SEE ALSO: MYFS_VIEWSURF  MYFS_VIEWSURFVOL
%
% (cc) 2020, dr.seunggoo.kim@gmail.com
if nargin < 2
  help(mfilename)
  return
end
if ~exist('cfg','var')
  cfg = [];
end

if ischar(data)
  % Source filename:
  fnameVolData = data;

  % Surface projection:
  [~, ~, surfData] = myfs_mni2fs(fnameVolData, struct(nofigure=1, fsavg='fsaverage5'));

  % Put them together:
  data = [surfData, fnameVolData];
end

if isempty(spaces)
  spaces = myfs_readsurfs('fsaverage5');
  spaces.mri = fullfile(getenv('SUBJECTS_DIR'),'cvs_avg35_inMNI152','mri','norm.mgz');
end

[H, cfg] = myfs_viewsurfvol(spaces, data, cfg);

end