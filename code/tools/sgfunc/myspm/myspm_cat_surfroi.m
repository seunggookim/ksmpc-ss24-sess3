function myspm_cat_surfroi(job)
% myspm_cat_checkvol(job)
%
% job
%  .dn_cat
%  .name_t1w
%
%

dn_cat = job.dn_cat;
name_t1w = job.name_t1w;
fn = ['lh.thickness.',name_t1w];
files = cellstr(spm_select('FPListRec',dn_cat, fn));

delete(gcp('nocreate')); parpool(8)
parfor i = 1:numel(files)
matlabbatch = {};
matlabbatch{1}.spm.tools.cat.stools.surf2roi.cdata = {files(i)};
spm_jobman('run',matlabbatch)
end
delete(gcp('nocreate'));
end