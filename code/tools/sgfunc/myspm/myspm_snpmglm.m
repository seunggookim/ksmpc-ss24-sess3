function JOB=myspm_snpmglm(JOB)
% check data
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
[~,~]=mkdir(JOB.dir_glm);

%%
matlabbatch={};
matlabbatch{1}.spm.tools.snpm.des.OneSampT.DesignName = ...
'MultiSub: One Sample T test on diffs/contrasts';
matlabbatch{1}.spm.tools.snpm.des.OneSampT.DesignFile = ...
'snpm_bch_ui_OneSampT';
matlabbatch{1}.spm.tools.snpm.des.OneSampT.dir = ...
{JOB.dir_glm};
matlabbatch{1}.spm.tools.snpm.des.OneSampT.P = fnames;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.cov = struct('c', {}, 'cname', {});
matlabbatch{1}.spm.tools.snpm.des.OneSampT.nPerm = JOB.nPerm;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.vFWHM = JOB.vFWHM;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.bVolm = 1;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.ST.ST_none = 0;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.masking.tm.tm_none = 1;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.masking.im = 1;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.masking.em = {''};
matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalc.g_omit = 1;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm.glonorm = 1;
matlabbatch{2}.spm.tools.snpm.cp.snpmcfg = ...
{[matlabbatch{1}.spm.tools.snpm.des.OneSampT.dir{1},'/SnPMcfg.mat']};

spm_jobman('initcfg')
spm_jobman('run', matlabbatch)
%% create "spmT_0001.nii" and "sigclus_+1.nii" and so on
cd (JOB.dir_glm)
alpha=0.05;
myunix(['fslmaths lP_FWE+.img -thr ',num2str(-log10(alpha)),' -bin sigclus_+1.nii'])
myunix(['fslmaths lP_FWE-.img -thr ',num2str(-log10(alpha)),' -bin sigclus_-1.nii'])
myunix(['fslmaths snpmT+.img spmT_0001.nii'])
myunix(['fslmaths snpmT-.img spmT_0002.nii'])
end
