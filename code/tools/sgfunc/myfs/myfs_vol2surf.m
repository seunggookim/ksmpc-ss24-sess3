function job = myfs_vol2surf (job)
% JOB = myfs_vol2surf (JOB) - This function is a wrapper of mri_vol2surf with 
% some aids for visual inspection of process results.
%
% mri_vol2surf samples the (registered) volume onto surfaces with a number of 
% sanity checks.
% If not registered, use fsss_bbr.m for inter-modality coregistation.
%
% JOB requires:
%  .fsdir      <1xN char> /path/where/you/have/freesurfer/subjects where you 
%                         have freesurfer subjects
%  .subjid     <1xN char> subject's ID in .fsdir
%  .subpath    <1xN char> where you want to save surface-mapped data
%  .fname_vol  <1xN char> input (functional/diffusion/stat...) volume
%  .regarg     <1xN char> filename of .dat registration file
%                         if the volume is in FS-space enter []
% (.projfrom)  <1xN str> 'white' (default)
% (.projarg)   <1xN str> e.g., " --projfrac 0.5" or " --projfrac-avg 0.1 0.9 0.1 "
% (.caxis)     <1x2 num> for coregistration
% (.interpopt)  'trilinear' (default) or 'nearest'
% (.outsuffix) <1xN str> output suffix
% (.fwhm_mm)   <1x1 num> surface-based smoothing (mm)
% (.optarg)
%
% JOB returns:
%  .
% (cc) 2015, 2017. sgKIM. solleo@gmail.com

if ~nargin, help (mfilename), return, end
if isfield(job,'overwrite'),overwrite=job.overwrite; else,overwrite=0;end
[~,myname] = fileparts(mfilename('fullpath'));
subjid = job.subjid;
disp(['@ Running ',myname,' on ',subjid,' with ',job.fname_vol,' ..']);
fsdir = job.fsdir;
setenv('SUBJECTS_DIR', fsdir);
subpath = job.subpath;
dir_exp = [fsdir,'/',subjid,'/',subpath];
if ~exist(dir_exp,'dir')
  [~,~] = mkdir(dir_exp);
end
if ~isfield(job,'projfrom'),  job.projfrom = 'white'; end
projfrom = job.projfrom;
if ~isfield(job,'interpopt')
  job.interpopt = 'trilinear';
end
interpopt = job.interpopt;
if ~isfield(job,'outsuffix'), job.outsuffix=subpath; end
if ~strcmp(job.outsuffix(1),'.')
  job.outsuffix = ['.',job.outsuffix];
end
outsuffix = job.outsuffix;
SIDE = {'lh','rh'};
[p1, f1, e1] = myfileparts(job.fname_vol);
fname_vol = [p1,'/',f1,e1];
ls(fname_vol);
if ~isfield(job,'regarg') || isempty(job.regarg)
  regarg = ['--regheader ',subjid];
else
  regarg = ['--reg ',job.regarg];
end
if ~isfield(job,'projarg'), job.projarg='--projfrac 0.5'; end
projarg = job.projarg;
Y = cell(1,2);
if ~isfield(job,'optarg'), job.optarg=''; end
if isfield(job,'fwhm_mm')
  outsuffix=[outsuffix,'.s',num2str(job.fwhm_mm)];
  job.optarg=[job.optarg,' --surf-fwhm ',num2str(job.fwhm_mm)];
end
for s=1:2
  fname_out=[dir_exp,'/',SIDE{s},outsuffix,'.mgz'];
  if ~exist(fname_out,'file') || overwrite
    unix(['mri_vol2surf --mov ',fname_vol,' ',regarg, ...
      ' --hemi ',SIDE{s},' --surf ',projfrom,' ',projarg, ...
      ' --sd ',fsdir,' --interp ',interpopt,' --o ', fname_out, ...
      ' ',job.optarg]);
  end
  ls(fname_out)
  Y{s} = squeeze(load_mgh(fname_out,[],1)); % for fMRI, #verts x 1 x 1 x #vols
  Y{s}(Y{s}==0) = nan;
  job.fname_out{s} = fname_out;
end

%% Surface visualization for sanity check
fname_png = [p1,'/',f1,'.',subjid,'.png'];
if ~isfield(job,'nofigure'), job.nofigure = 0; end
if (~exist(fname_png,'file') || overwrite) && (~job.nofigure)
  surfs = fsss_read_all_FS_surfs(subjid, fsdir);
  cfg = struct('colorbartitle',[subjid,outsuffix], 'fname_png',fname_png,...
    'dpi',200, 'basesurf','white', 'colorbarinterp','none',...
    'colormap',parula);
  if isfield(job,'caxis'), cfg.caxis = job.caxis; end
  if isfield(job,'colormap'), cfg.colormap = job.colormap; end
  myfs_view(surfs, Y, cfg);
end
disp(['> Done:',myname,' on ',subjid])
end
