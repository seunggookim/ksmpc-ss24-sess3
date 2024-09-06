function mri = helper_readmri(fname)

[~,~,ext] = fileparts_gz(fname);
switch (ext)
  case {'.nii','.nii.gz','.img'}
    % Try "my" version of MATLAB Image Processing Toolbox (since 2016)
    if exist('niftiinfogz','file') && exist('niftireadgz','file')
      [V,info] = niftireadgz(fname);
      mri = struct('vol',V, 'info',info);
      % NOTE: NIFTIREAD doesn't upscale precision (hmm...)
    elseif exist('load_nifti','file') % FreeSurfer
      mri = load_nifti(fname);
    elseif exist('load_untouch_nii','file') % NIFTI toolbox
      mri = load_untouch_nii(fname);
    else
      error('CANNOT FIND any function to read NIFTI/ANALYZE files!')
    end

  case {'.mgh','.mgz'}
    [~, M, P] = load_mgh(fname, [], [], 1);
    mri = struct('vox2ras',M, 'tr',P(1), 'filpangle',P(2), 'te',P(3), 'ti',P(3), 'fov',P(4));
    [mri.vol] = load_mgh(fname);

  otherwise
    error('EXTENSION UNRECOGNIZED: %s',ext)
end

mri.vol = double(mri.vol); % now I just upscale it as I won't save it
end