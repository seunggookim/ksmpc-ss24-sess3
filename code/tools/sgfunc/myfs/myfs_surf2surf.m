function [job,Y,Cmd] = myfs_surf2surf(job)
% [JOB,Y] = myfs_surf2surf(JOB)
%
% JOB requires:
%  .fname_suffix [1xN] lh.${fname_suffix} 
%  .fsdir        '1xN'
%  .subpath      '1xN'
%  .srcsubj      '1xN'
%  .trgsubj      '1xN'
% (.fwhm_mm)     [1x1]
% (.mapmethod)   '1xN' nnfr (default: linear) | nnf (nearest)
%
%
% (cc) 2015,2019, sgKIM

if isfield(job,'overwrite'),overwrite=job.overwrite;else,overwrite=0;end
[~,myname] = fileparts(mfilename('fullpath'));
SIDE = {'lh','rh'};
fsdir = job.fsdir;
subpath = job.subpath;
setenv('SUBJECTS_DIR',fsdir);
srcsubj = job.srcsubj;
trgsubj = job.trgsubj;
disp(['@ Running ',myname,': ',srcsubj,' to ',trgsubj,' ..']);
[~,~] = mkdir([job.fsdir,'/',trgsubj,'/',subpath]);
fwhmarg = ''; 
fwhmsuffix = '';
if isfield(job,'fwhm_mm')
  fwhmarg = [' --fwhm-trg ',num2str(job.fwhm_mm),' '];
  fwhmsuffix = ['.s',num2str(job.fwhm_mm)];
end
[~,f1,e1] = fileparts(job.fname_suffix);
suffix1 = job.fname_suffix;
if contains(e1, {'mgh','mgz','w'})
  suffix2 = [f1,'.',trgsubj];
else
  suffix2 = [f1,e1,'.',trgsubj];
end
if isfield(job, 'mapmethod')
  mapmethod = [' --mapmethod ',job.mapmethod];
else
  mapmethod = '';
end
Cmd = {};

for s = 1:2
  prefix = [job.fsdir,'/',srcsubj,'/',subpath,'/', SIDE{s}];  
  input = [prefix,'.',suffix1];
  output = [prefix,'.',suffix2,fwhmsuffix,'.mgz'];
  cmd = [getenv('FREESURFER_HOME'),'/bin/mri_surf2surf ', ...
    ' --srcsubject ',srcsubj,' --sval ',input, ...
    ' --trgsubject ',trgsubj,' --tval ',output, ...
    ' --hemi ',SIDE{s}, fwhmarg,' ',mapmethod];
  Cmd = [Cmd cmd];
end
Y = {};  
if isfield(job,'GenCmdOnly') && job.GenCmdOnly
  return
end

for s = 1:2
  prefix = [job.fsdir,'/',srcsubj,'/',subpath,'/', SIDE{s}];  
  input = [prefix,'.',suffix1];
  output = [prefix,'.',suffix2,fwhmsuffix,'.mgz'];
  if ~isfile(output) || overwrite
    system(Cmd{s})
  end
  ls(output)
  Y{s}=squeeze(load_mgh(output,[],1));
end

%% Surface visualization for sanity check
dname_fig = [fsdir,'/',srcsubj,'/',subpath,'/fig/'];
fname_png=[dname_fig,suffix2,fwhmsuffix,'.png'];
if ~isfield(job,'nofigure'), job.nofigure = 0; end
if (~isfile(fname_png) || overwrite) && (~job.nofigure)
  [~,~]=mkdir(dname_fig);
  surfs = myfs_readsurfs(trgsubj, fsdir);
  cfg = struct('colorbartitle',[srcsubj,'.',suffix2,fwhmsuffix],...
    'fname_png',fname_png, 'dpi',100, 'basesurf','white', ...
    'colorbarinterp','none', 'colormap',parula);
  if isfield(job,'caxis'), cfg.caxis = job.caxis; end
  if isfield(job,'colormap'), cfg.colormap = job.colormap; end
  map = {Y{1}(:,1), Y{2}(:,1)};
  map{1}(map{1}==0) = nan;
  map{2}(map{2}==0) = nan;
  myfs_view(surfs, map, cfg);
end
disp(['> Done:',myname,': ',srcsubj,' to ',trgsubj]);
end
