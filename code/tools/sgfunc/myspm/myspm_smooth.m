function JOB = myspm_smooth(JOB)
% JOB = myspm_smooth(JOB)
%
% JOB requires:
%  .fname    '1xN' for a file name
%  .fwhm_mm  [1x1] or [1x3]
%
% (cc) 2015, sgKIM.

if ~nargin, help(mfilename); return; end
if ~isfield(JOB,'overwrite'), JOB.overwrite = 0; end

ls(JOB.fname);
matlabbatchs = {};
matlabbatchs{1}.spm.spatial.smooth.dtype = 0;  % same precision
matlabbatchs{1}.spm.spatial.smooth.im = 0;     % implicit masking
%% Find # of volumes and set them as inputs
hdr = spm_vol_nifti(JOB.fname);
if numel(hdr.private.dat.dim) > 3
  T = hdr.private.dat.dim(4);
else
  T = 1;
end
matlabbatchs{1}.spm.spatial.smooth.data = cell(T,1);
for t = 1:T
  matlabbatchs{1}.spm.spatial.smooth.data{t,1} ...
    = [JOB.fname,',',num2str(t)];
end
%% Set a prefix for the kernel size
if numel(JOB.fwhm_mm) == 3 % input is [1x3] or [3x1]
  if (JOB.fwhm_mm(1) == JOB.fwhm_mm(2)) ...
      && (JOB.fwhm_mm(2) == JOB.fwhm_mm(3)) % isotropic?
    f = num2str(JOB.fwhm_mm(3));
  else % anisotropic
    f = [num2str(JOB.fwhm_mm(1)),'x',num2str(JOB.fwhm_mm(2)), ...
      'x',num2str(JOB.fwhm_mm(3))];
  end
elseif numel(JOB.fwhm_mm)==1 % input is [1x1]
  f=num2str(JOB.fwhm_mm);
  JOB.fwhm_mm = [1 1 1].*JOB.fwhm_mm;
else
  error(['# JOB.fwhm_mm must be either a scalar or [1x3] vector!']);
end
matlabbatchs{1}.spm.spatial.smooth.fwhm = JOB.fwhm_mm;
prefix = ['s',f];
matlabbatchs{1}.spm.spatial.smooth.prefix = prefix;
%% run if needed
[p1,f1,e1] = myfileparts(JOB.fname);
fname_out = fullfile(p1,['/' prefix f1 e1]);
if ~exist(fname_out,'file') || JOB.overwrite
  spm_jobman('initcfg')
  spm_jobman('run', matlabbatchs)
end
end
