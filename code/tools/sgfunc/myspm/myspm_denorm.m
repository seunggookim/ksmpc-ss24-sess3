function JOB = myspm_denorm(JOB)
% JOB = myspm_denorm(JOB)
%
% applying inversed normalization (putting back to the native space)
% saving the outputs in the source directory (can be changed)
%
% JOB requires:
%  .fname_fixed
%  .fname_moving
%  .fname_deform (y_*)
% (.interp) B-spline order from 0 to 6 (default=4)
% (.prefix) default='i'
% (.savepwd)
%
% (cc) 2018, sgKIM, solleo@gmail.com

if nargin==0, help(mfilename); return; end
spm('Defaults','fmri');
% this is SPM12-included batch process using SPM12
a=spm('version');
if ~strcmp(a(4:5),'12'),  error(['Run ',mfilename,' on SPM12!']);  end
if ~isfield(JOB,'interp'), JOB.interp = 4; end % 4-spline as default

matlabbatch={};
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {JOB.fname_deform};
matlabbatch{1}.spm.util.defs.comp{1}.inv.space  = {JOB.fname_fixed};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {JOB.fname_moving};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = JOB.interp;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
if isfield(JOB,'savepwd') && JOB.savepwd
  matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savepwd = 1;
else
  matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
end
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'i';
%save(fullfile(pwd,['spm12_denorm.mat']), 'matlabbatch');
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);

end
