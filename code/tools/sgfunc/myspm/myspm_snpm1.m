function myspm_snpm1(JOB)

% initialize SPM
spm('Defaults','fmri')

% find data
if isfield(JOB,'fnames')
fnames = JOB.fnames;
elseif isfield(JOB,'files_query')
JOB.subjID=fsss_subjID(JOB.subjID);
fnames=cell(numel(JOB.subjID),1); % this MUST be a column vector <Nx1>
idx = strfind(JOB.files_query,'${subj}');
prefix = JOB.files_query(1:idx-1);
suffix = JOB.files_query(idx+7:end);
for n=1:numel(JOB.subjID)
[~,res] = mydir([prefix,JOB.subjID{n},suffix]);
if isempty(res)
error(['File not found: ',prefix,JOB.subjID{n},suffix]);
end
fnames{n,1}=[res,',1'];
end
elseif isfield(JOB,'filenames')
for n=1:size(JOB.filenames,1)
fnames{n,1} = [JOB.filenames{n,1},',1'];
end
else
error('You need to specify inputs in JOB.files_query or JOB.filenames');
end
Nsubj=numel(fnames);
for n=1:Nsubj
ls(fnames{n,1}(1:end-2));
end

% find dir
if ~isfield(JOB,'dir_glm')%&&isfield(JOB,'dir_base');
if ~isfield(JOB,'dir_prefix'), JOB.dir_prefix=''; end
JOB.dir_glm=fullfile(JOB.dir_base,[JOB.dir_prefix,JOB.model_desc]);
end
[~,~]=mkdir(JOB.dir_glm);

% set spnm1
snpm1.DesignName = 'MultiSub: One Sample T test on diffs/contrasts';
snpm1.DesignFile = 'snpm_bch_ui_OneSampT';
snpm1.dir = {JOB.dir_glm};
snpm1.P = fnames;
snpm1.cov = struct('c', {}, 'cname', {});
snpm1.nPerm = 5000;
snpm1.vFWHM = [0 0 0];
snpm1.bVolm = 1;
snpm1.ST.ST_none = 0;
snpm1.masking.tm.tm_none = 1;
snpm1.masking.im = 1;
snpm1.masking.em = {''};
snpm1.globalc.g_omit = 1;
snpm1.globalm.gmsca.gmsca_no = 1;
snpm1.globalm.glonorm = 1;

matlabbatch={};
matlabbatch{1}.spm.tools.snpm.des.OneSampT=snpm1;

save /tmp/mat.mat matlabbatch

spm_jobman('initcfg')
spm_jobman('run', matlabbatch)
%cd(JOB.dir_glm)
%spm_print;


end
