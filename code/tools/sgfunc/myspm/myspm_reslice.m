function JOB = myspm_reslice(JOB)
% JOB requires
%  .fname_ref
%  .fname_src
% (.interp)
% (.prefix)
%
% (cc) 2017, sgKIM
if ~nargin, help(mfilename); return; end
if isfield(JOB,'interp')
 interp=JOB.interp;
else
 interp=4;
end
matlabbatch={};
matlabbatch{1}.spm.spatial.coreg.write.ref = {[JOB.fname_ref,',1']};
matlabbatch{1}.spm.spatial.coreg.write.source = ...
  cellstr(spm_select('expand',JOB.fname_src));
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = interp;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
if ~isfield(JOB,'prefix')
  JOB.prefix='r';
end
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = JOB.prefix;

spm_jobman('initcfg')
spm_jobman('run', matlabbatch)
[p1,f1,e1] = myfileparts(JOB.fname_src);
JOB.fname_out = [p1,'/',JOB.prefix,f1,e1];

end