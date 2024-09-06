function [fname_out, surfs, Y] = myfs_mni2fs(fname_in, cfg)
% [fname_out, surfs, Y] = myfs_mni2fs(fname_in, cfg)
%
%
% fname_in : NIFTI file in MNI152 space
%
% cfg options:
%   .fsavg      Freesurfer template surface (Default=fsaverage6)
%   .interp     'nearest' (default)
%
% (cc) 2021, sgKIM.

%{
NB: mri_vol2surf can project to trgsubject but nearest interpolation is
not working. So here I use mri_vol2surf and mri_surf2surf consecutively
(it doesn't take not too much longer).
%}

[p1,f1,~] = myfileparts(fname_in);
if ~exist('cfg','var'), cfg = []; end
if ~isfield(cfg,'fsavg')
  cfg.fsavg = 'fsaverage6';
end
if ~isfield(cfg,'interp')
  cfg.interp = 'nearest';
end
switch cfg.interp
  case 'nearest'
    cfg.interopt = 'nearest';
    cfg.mapmethod = 'nnf';
  case 'trilinear'
    cfg.interopt = 'trilinear';
    cfg.mapmethod = 'nnfr';
end
fname_out = [p1,'/',f1,'.',cfg.fsavg,'.mat'];
if isfield(cfg,'suffix')
  fname_out = [p1,'/',f1,'.',cfg.fsavg,'.',cfg.suffix,'.mat'];
end
srcsubj = 'cvs_avg35_inMNI152';
hemis = {'lh','rh'};
if ~exist(fname_out,'file')
  for i = 1:2
    fname_surf{i} = [tempname,'.mgh'];
    cmd = sprintf(...
      'mri_vol2surf --mov %s --hemi %s --regheader %s --o %s --interp nearest',...
      fname_in, hemis{i}, srcsubj, fname_surf{i});
    unix(cmd)
    cmd = sprintf(...
      'mri_surf2surf --srcsubject %s --trgsubject %s --sval %s --tval %s --mapmethod %s --hemi %s',...
      srcsubj, cfg.fsavg, fname_surf{i}, fname_surf{i}, cfg.mapmethod, hemis{i});
    unix(cmd)
  end
  Y = cell(1,2);
  for i = 1:2
    Y{i} = squeeze(load_mgh(fname_surf{i},[],1));
  end
  save(fname_out,'Y');
end
if (~isfield(cfg,'nofigure') || ~cfg.nofigure) || (nargout>1)
  surfs = myfs_readsurfs(cfg.fsavg);
  if ~exist('Y','var'), load(fname_out,'Y'); end
end
if (~isfield(cfg,'nofigure') || ~cfg.nofigure)
  myfs_view(surfs, Y, cfg);
end

end
