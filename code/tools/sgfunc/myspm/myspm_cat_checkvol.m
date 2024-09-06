function myspm_cat_checkvol(job)
% myspm_cat_checkvol(job)
%
% job
%  .dn_cat
%  .name_t1w
%
%

% dn_cat = '/home/sgk/github/HL-t1w/data/cat12/';
% name_t1w = 't1w';
dn_cat = job.dn_cat;
name_t1w = job.name_t1w;
dn_fig = fullfile(dn_cat,'fig');
mkdir (dn_fig)

%% First view one slice from ALL subjects:
fn = ['wm',name_t1w,'.nii'];
files = cellstr(spm_select('ExtFPListRec',dn_cat, fn, 1));
matlabbatch = {};
matlabbatch{1}.spm.tools.cat.tools.showslice.data_vol = files;  
matlabbatch{1}.spm.tools.cat.tools.showslice.scale = 0;
matlabbatch{1}.spm.tools.cat.tools.showslice.orient = 3;
matlabbatch{1}.spm.tools.cat.tools.showslice.slice = 0;
spm_jobman('run',matlabbatch)

spm_print(fullfile(dn_fig,['all_wm',name_t1w,'.ps']));

%% [INTERACTIVE] Sample homogeniety: volume
% ext-fullpath-list-recursive (for images):
fn = ['wm',name_t1w,'.nii'];
niis = cellstr(spm_select('ExtFPListRec', dn_cat, fn, 1));
% fullpath-list-recursive (for other files):
fn = ['cat_',name_t1w,'.xml'];
xmls = cellstr(spm_select('FPListRec', dn_cat, fn));

matlabbatch = {};
matlabbatch{1}.spm.tools.cat.tools.check_cov.data_vol = {niis}; % group1
matlabbatch{1}.spm.tools.cat.tools.check_cov.data_xml = xmls;
matlabbatch{1}.spm.tools.cat.tools.check_cov.gap = 3;
matlabbatch{1}.spm.tools.cat.tools.check_cov.c = {};
spm_jobman('run',matlabbatch)
spm_print(fullfile(dn_fig,['vol_meancorr.ps']));

end