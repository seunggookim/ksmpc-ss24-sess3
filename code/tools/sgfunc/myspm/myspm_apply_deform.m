function myspm_apply_deform(JOB)
% JOB = myspm_apply_deform(JOB)
%
% JOB requires:
% .fname_deform
% .isInv
% .fname_moving
% .fname_fixed
% .intorder
%
% (cc) 2017, sgKIM.

if ~isfield(JOB,'intorder'), intorder=1; else, intorder=JOB.intorder; end
if ~isfield(JOB,'isInv'), isInv=0; else, isInv=JOB.isInv; end

matlabbatch={};
if isInv
 matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {JOB.fname_deform};
 matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {JOB.fname_fixed};
else
 matlabbatch{1}.spm.util.defs.comp{1}.def = {JOB.fname_deform};
end
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {JOB.fname_moving};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savepwd = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = intorder;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask   = 0;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm   = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.prefix = 'd';

spm_jobman('initcfg');
spm_jobman('run', matlabbatch);
end
