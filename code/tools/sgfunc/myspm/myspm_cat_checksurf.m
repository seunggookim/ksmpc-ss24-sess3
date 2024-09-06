function myspm_cat_checksurf(job)
% myspm_cat_checkvol(job)
%
% job
%  .dn_cat
%  .name_t1w
%
%
dn_fig = fullfile(job.dn_cat,'fig');
[~,~] = mkdir (dn_fig);

%% Check if thickness is resampled:

fn = ['^',job.name_t1w,'.nii'];
files = cellstr(spm_select('FPListRec', job.dn_cat, fn));
nsub = numel(files);

fn = ['mesh.thickness.resampled_32k.',job.name_t1w,'.gii'];
files = cellstr(spm_select('FPListRec', job.dn_cat, fn));
if numel(files) ~= nsub
  myspm_cat_surf_resample(job)
end

%% [INTERACTIVE] Sample homogeniety: thickness

% fullpath-list-recursive (for other files):
fn = ['mesh.thickness.resampled_32k.',job.name_t1w,'.gii'];
files = cellstr(spm_select('FPListRec', job.dn_cat, fn));
% fullpath-list-recursive (for other files):
fn = ['cat_',job.name_t1w,'.xml'];
xmls = cellstr(spm_select('FPListRec', job.dn_cat, fn));

matlabbatch = {};
matlabbatch{1}.spm.tools.cat.stools.check_mesh_cov.data_surf = {files}; % group1
matlabbatch{1}.spm.tools.cat.stools.check_mesh_cov.data_xml = xmls;
matlabbatch{1}.spm.tools.cat.stools.check_mesh_cov.c = {};

spm_jobman('run',matlabbatch)
spm_print(fullfile(dn_fig,['surf_meancorr.ps']));

end