function myspm_cat_surf_resample(job)
% myspm_cat_surf_resample(job)
%
% job
%  .dn_cat
%  .name_t1w
% (.name_metric = 'thickness')
% (.fwhm_surf   = 0)

% merge_hemi?
% mesh32k?
% nproc?

if ~isfield(job,'name_metric')
  job.name_metric = 'thickness';
end
if ~isfield(job,'fwhm_surf')
  job.fwhm_surf = 0;
end
if ~isfield(job,'nproc')
  job.nproc = 8;
end

fn = ['lh.',job.name_metric,'.',job.name_t1w];
files = cellstr(spm_select('FPListRec', job.dn_cat, fn));
delete(gcp('nocreate')); parpool(job.nproc)
parfor i = 1:numel(files)
  matlabbatch = {};
  matlabbatch{1}.spm.tools.cat.stools.surfresamp.data_surf = files(i);
  matlabbatch{1}.spm.tools.cat.stools.surfresamp.merge_hemi = 1;
  matlabbatch{1}.spm.tools.cat.stools.surfresamp.mesh32k = 1;
  matlabbatch{1}.spm.tools.cat.stools.surfresamp.fwhm_surf = job.fwhm_surf;
  matlabbatch{1}.spm.tools.cat.stools.surfresamp.nproc = 0; % now sure why but this won't run pretty: it's trying to do N^2 times...
  
  [dn,~,~] = fileparts(files{i});
  % well... for now... (only for mesh32k=1, merge_hemi=1)
  fn = ['mesh.',job.name_metric,'.resampled_32k.',job.name_t1w,'.gii'];
  if job.fwhm_surf
    fn = ['s',num2str(job.fwhm_surf) fn];
  end
  fn_out = fullfile(dn, fn);
  
  if ~isfile(fn_out)
    spm_jobman('run',matlabbatch)
  end
end
delete(gcp('nocreate'));
end