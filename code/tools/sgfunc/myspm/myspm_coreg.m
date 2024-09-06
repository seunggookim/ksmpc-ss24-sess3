function JOB = myspm_coreg(JOB)
% JOB = myspm_coreg(JOB)
%
% transforms source (moving/others) files without modifying headers
%
% JOB requires:
%  .prefix
%  .interp
%  .fname_moving
%  .fname_fixed
%  .fname_others
%
% (cc) 2015, 2019, sgKIM.
if ~nargin,  help myspm_coreg;  return; end
if ~isfield(JOB,'interp'), JOB.interp=4; end
if ~isfield(JOB,'overwrite'); JOB.overwrite=0; end
if ~isfield(JOB,'prefix'), JOB.prefix='o'; end

[p1,f1,e1] = myfileparts(JOB.fname_moving);
JOB.fname_moving = [p1,'/',f1,e1];
[p2,f2,e2] = myfileparts(JOB.fname_fixed);
JOB.fname_fixed = [p2,'/',f2,e2];

ls(JOB.fname_moving);
ls(JOB.fname_fixed);
% make a temp source file:
fname_moving_temp = [tempname,e1];
copyfile(JOB.fname_moving, fname_moving_temp, 'f');
fname_moving{1,1} = [fname_moving_temp,',1'];
pwd0=pwd;
cd(p1);
fnames_others={};
if isfield(JOB,'fname_others')
  if ~iscell(JOB.fname_others)
    JOB.fname_others={JOB.fname_others};
  end
  for c=1:numel(JOB.fname_others)
    ls(JOB.fname_others{c});
    % make temp source files:
    [p3,f3,e3] = myfileparts(JOB.fname_others{c});
    fname_others_temp{c} = [tempname,e3];
    copyfile(JOB.fname_others{c}, fname_others_temp{c}, 'f');
    fnames_others=[fnames_others; ...
      cellstr(spm_select('expand',fname_others_temp{c}))];
  end
end

matlabbatch = {};
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[JOB.fname_fixed,',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = fname_moving;
if isfield(JOB,'fname_others')
  matlabbatch{1}.spm.spatial.coreg.estwrite.other = fnames_others;
end
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol ...
  = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = JOB.interp;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = JOB.prefix;
save([p1,'/coreg_',datestr(now,'yyyymmdd'),'.mat'], 'matlabbatch');
spm('Defaults','fmri')
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);

% delete temp files
delete(fname_moving_temp)
if isfield(JOB,'fname_others')
  for c=1:numel(JOB.fname_others)
    delete(fname_others_temp{c})
  end
end

% change output filenames
[p4,f4,e4] = myfileparts(fname_moving_temp);
movefile([p4,'/',JOB.prefix,f4,e4], [p1,'/',JOB.prefix,f1,e1])
if isfield(JOB,'fname_others')
  for c=1:numel(JOB.fname_others)
    [p4,f4,e4] = myfileparts(fname_others_temp{c});
    [p3,f3,e3] = myfileparts(JOB.fname_others{c});
    movefile([p4,'/',JOB.prefix,f4,e4], [p3,'/',JOB.prefix,f3,e3])
  end
end

cd(pwd0)
