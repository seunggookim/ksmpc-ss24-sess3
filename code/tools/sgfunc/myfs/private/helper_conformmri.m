function mri = helper_conformmri(mri)
if ~isfield(mri,'vol')
  if isfield(mri,'img')  % from LOAD_UNTOUCH_NII
    mri.vol = mri.img;
    mri = rmfield(mri,'img');
  else
    error('UNKNOWN volume fieldname??')
  end
end
if ~isfield(mri,'vox2ras')
  if isfield(mri,'hdr')  % from LOAD_UNTOUCH_NII
    mri.vox2ras = [mri.hdr.hist.srow_x; mri.hdr.hist.srow_y; mri.hdr.hist.srow_z; 0 0 0 1];
  elseif isfield(mri,'info')  % from NIFTIINFO
    mri.vox2ras = mri.info.Transform.T';  % 0-based
  else
    error('UNKNOWN header fieldname??')
  end
end
if size(mri.vol,4) > 1
  fprintf('[%s] 4-D image is given. Showing the only first volume.\n', ...
    mfilename);
  mri.vol = mri.vol(:,:,:,1);
end
mri.vol = double(mri.vol);
end