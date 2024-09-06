function myspm_unnorm(JOB)
% myspm_unnorm(JOB)
%
% JOB requires:
%  .fname_y        : T1w deformation field (e.g. y_t1w.nii)
%  .fname_ref      : EPI coregistered in native-T1w space (e.g., uaEPI.nii)
%  .fname_input    : Image in MNI152 space
% (.prefix)        : (default=x)
% (.interp)        : (default=4; 4-th order B-spline)

defs=[];
ls(JOB.fname_y)
defs.comp{1}.inv.comp{1}.def = {JOB.fname_y};
ls(JOB.fname_ref)
defs.comp{1}.inv.space = {JOB.fname_ref};
ls(JOB.fname_input)
hdr=load_untouch_header_only(JOB.fname_input);
% T=hdr.dime.dim(5);
% for t=1:T
% defs.out{1}.pull.fnames{t,1} = [JOB.fname_input,',',num2str(t)];
% end
defs.out{1}.pull.fnames{1}=JOB.fname_input;
defs.out{1}.pull.savedir.savesrc = 1;
if ~isfield(JOB,'interp'), JOB.interp=4; end
defs.out{1}.pull.interp = JOB.interp;
defs.out{1}.pull.mask = 1;
defs.out{1}.pull.fwhm = [0 0 0];
if ~isfield(JOB,'prefix'), JOB.prefix='x'; end
defs.out{1}.pull.prefix = JOB.prefix;
matlabbatch={};
matlabbatch{1}.spm.util.defs = defs;
spm_jobman('run',matlabbatch)
[p1,f1,e1]=fileparts(JOB.fname_input);
fname_out=[p1,'/',JOB.prefix,f1,e1];
ls(fname_out)
myspm_denan(fname_out);
end