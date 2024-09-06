function job = myfs_smooth(job)
%MYFS_SMOOTH smoothes data within the cortex/brain (fs-fast).
%
% [USAGE]
% job = MYFS_READSURFS (job)
%
% [INPUT]
% job (1x1) job description structure:
% - space:
%  .subid        '1xN' subject directory name
%  .subjects_dir '1xN' directory where you have the subject directory
%
% - input: (two possible options)
% -- A. MGZ files
%  .fname_src    '1xN' any surface mapped data in .MGZ ("?h.xxx.mgz")
% -- B. FS-FAST
%  .fsd          '1xN' assumes the data structure processed FS-FAST
%  .fstemp       '1xN' template surfaces that the data is normalized into
% (.mni305)      [1x1] 1=1mm | 2=2mm (default)
%  .runs         [1xN] runs to smooth
%
% - smoothing:
%  .fwhm_mm      [1x1] full-width at half maximum in mm
%
% - output:
% (.overwrite)   [1x1] true=overwrite | false=don't (default)
%
% (cc) 2020, sgKIM, solleo@gmail.com


%% CHECK inputs & SET default values
if ~nargin
  help(mfilename)
  return
end

if isfield(job,'subjects_dir')
  setenv('SUBJECTS_DIR',job.subjects_dir);
else
  job.subjects_dir = getenv('SUBJECTS_DIR');
end
if ~isfield(job,',mni305')
  job.mni305 = 2;
end
if ~isfield(job,'overwrite')
  job.overwrite = false;
end

%% Running

if isfield(job,'fname_src') && ~isfield(job,'fsd')
  % INPUT A: MGZ FILES
  dname = fullfile(job.subjects_dir, job.subid);
  hemis = {'lh','rh'};
  for ihemi = 1:2
    hemi = hemis{ihemi};
    [path1,fname1,ext1] = fileparts_gz(job.fname_src);
    if isempty(path1)  % if relative path, assuming it's in <surf>
      path1 = fullfile(dname,'surf');
    end
    src = fullfile(path1,[hemi,fname1(3:end),ext1]);
    trg = fullfile(path1,sprintf('%s%s.s%g%s',...
      hemi,fname1(3:end),job.fwhm_mm,ext1));
    if ~isfile(trg) || job.overwrite
      system(sprintf(...
        ['mri_surf2surf --cortex --hemi %s --s %s --sval %s --tval %s',...
        ' --fwhm-trg %g'], hemi, job.subid, src, trg, job.fwhm_mm));
    else
      ls(trg)
    end
  end
  
elseif ~isfield(job,'fname_src') && isfield(job,'fsd')
  % INPUT B: FS-FAST FILES
  for runid = job.runs
    dname = fullfile(job.subjects_dir, job.subid, job.fsd, ...
      sprintf('%03i',runid));
    
    % - Cortex
    hemis = {'lh','rh'};
    for ihemi = 1:2
      hemi = hemis{ihemi};
      src = fullfile(dname, sprintf('fmcpr.sm0.%s.%s.nii.gz', ...
        job.fstemp, hemi));
      trg = fullfile(dname, sprintf('fmcpr.sm%g.%s.%s.nii.gz',...
        job.fwhm_mm, job.fstemp, hemi));
      if ~isfile(trg) || job.overwrite
        system(sprintf(...
          ['mri_surf2surf --cortex --hemi %s --s %s --sval %s --tval %s',...
          ' --fwhm-trg %g'], hemi, job.fstemp, src, trg, job.fwhm_mm));
      else
        ls(trg)
      end
    end
    
    % - Subcortical
    sctx = sprintf('mni305.%imm',job.mni305);
    src = fullfile(dname, sprintf('fmcpr.sm0.%s.nii.gz', sctx));
    trg = fullfile(dname, sprintf('fmcpr.sm%g.%s.nii.gz',...
      job.fwhm_mm, sctx));
    mask = fullfile(dname, 'masks', sprintf('brain.%s.nii.gz', sctx));
    if ~isfile(trg) || job.overwrite
      system(sprintf(...
        ['mri_fwhm --i %s --o %s --smooth-only --fwhm %g --mask %s'],...
        src, trg, job.fwhm_mm, mask));
    else
      ls(trg)
    end
  end
  
else
  error('ENTER EITHER .job_fname_src (MGZ) or .fsd (FSFAST)!')
end





end

function [path,name,ext]=fileparts_gz(filename)
% [path,name,ext]=fileparts_gz(filename)
% (cc) 2014, sgKIM.

[path,name,ext] = fileparts(filename);
if strcmp(ext,'.gz')
  ext = [name(end-3:end),ext];
  name = name(1:end-4);
end

end
