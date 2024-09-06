function JOB = myspm_coreg4d(JOB)
% JOB = myspm_coreg(JOB)
%
% JOB requires:
%  .prefix
%  .interp
%  .fname_moving
%  .fname_fixed
%
% (cc) 2015, sgKIM.
if ~nargin,  help myspm_coreg;  return; end
if ~isfield(JOB,'interp'), JOB.interp=4; end
if ~isfield(JOB,'overwrite'); JOB.overwrite=0; end
if ~isfield(JOB,'prefix'), JOB.prefix='o'; end
fnames_moving{1,1} = [JOB.fname_moving,',1'];
ls(JOB.fname_moving);
ls(JOB.fname_fixed);
[p1,f1,e1]=fileparts(JOB.fname_moving);
if isempty(p1), p1=pwd; end
cd(p1);
% hdr = load_untouch_header_only(JOB.fname_moving);
% NT = hdr.dime.dim(5);
% fnames_others = cell(NT,1);
% for t=1:NT
%  fnames_others{t,1}=[JOB.fname_moving,',',num2str(t)];
% end
fnames_others=cellstr(spm_select('expand',JOB.fname_moving));
fname_out=fullfile(p1,[JOB.prefix,f1,e1]);
if ~exist(fname_out,'file') || JOB.overwrite
 matlabbatch={};
 matlabbatch{1}.spm.spatial.coreg.estwrite.ref    = {[JOB.fname_fixed,',1']};
 matlabbatch{1}.spm.spatial.coreg.estwrite.source = fnames_moving;
 matlabbatch{1}.spm.spatial.coreg.estwrite.other  = fnames_others;
 matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
 matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
 matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol ...
  = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
 matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm   = [7 7];
 matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = JOB.interp;
 matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap   = [0 0 0];
 matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask   = 0;
 matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = JOB.prefix;
 save([p1,'/coreg_',datestr(now,'yyyymmdd'),'.mat'], 'matlabbatch');
 spm('Defaults','fmri')
 spm_jobman('initcfg');
 spm_jobman('run', matlabbatch);
end
end
