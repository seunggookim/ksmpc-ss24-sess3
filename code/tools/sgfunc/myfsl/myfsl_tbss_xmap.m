function JOB = myfsl_tbss_xmap(JOB)
% JOB = myfsl_tbss_xmap(JOB)
%
% JOB requires:
%  .subjID
%  .dir_input  : all inputs should be ${dir_input}/${subj}${suffix}.nii
%  .dir_fa
%  .dir_output
%  .suffix...?
%
% NOTE:
% for some reasons... tbss_non_FA doesn't work well with probabilistic
% tractogram (vmap, nmap, cmap,... whatever)

if ~nargin,   help myfsl_tbss_xmap;  return;  end
if ~isfield(JOB,'alpha'),   JOB.alpha=0.05;  end
subjID= fsss_subjID(JOB.subjID);
N = numel(subjID);
if ~isfield(JOB,'glm_prefix'),  JOB.glm_prefix=''; end
setenv('FSLOUTPUTTYPE','NIFTI_GZ');
if isfield(JOB,'dir_base')
dir_base=JOB.dir_base;
else
dir_base=pwd;  
end

%% coreg (only because this is how I do in myfsl_tbss
dir1=[JOB.dir_base,'/inputs/'];
[~,~]=mkdir(dir1);
if ~isfield(JOB,'MEAS'),  MEAS={'cmap'}; else MEAS=JOB.MEAS; end
M=numel(MEAS);
for j=1:M
meas = MEAS{j};
dir3=[JOB.dir_base,'/',meas,'/'];
[~,~]=mkdir(dir3);
% visual inspection! (am i putting in right things?)
Ori={'sag','cor','axl'};
dir_fig = [JOB.dir_base,'/fig_overview/',meas,'/'];
[~,~]=mkdir(dir_fig);
if ~exist([dir_fig,'probtrck_',meas,'_00_native_',Ori{end},'.png'],'file')
job1=[];
job1.subjID = subjID;
job1.base_query=[dir1,'${subj}/t1w.nii'];
job1.over_query=[dir1,'${subj}/',meas,'.nii'];
job1.seed_query=[dir1,'${subj}/seed.nii'];
for k=1:3
job1.cfg=struct('slicedim',k, 'thres',0);
job1.fname_fig=[dir_fig,'probtrck_',meas,'_00_native_',Ori{k},'.png'];
imagethumbnail (job1);
end
end

IMG={'seed','t1w',meas};
M=numel(IMG);
for n=1:N
subjid=subjID{n};
job1.name_moving=[dir1,subjid,'/FA.nii'];
job1.name_fixed='/home/raid2/skim/FMRIB58_FA_1mm.nii';
job1.name_others={};
for j=1:M
if ~exist([dir1,subjid,'/r',IMG{j},'.nii'],'file')
job1.name_others={ job1.name_others{:},[dir1,subjid,'/',IMG{j},'.nii'] };
end
end
job1.interp=1;
job1.prefix='r';
if ~isempty(job1.name_others)
myspm_coreg(job1);
end
end

% visual inspection!
if ~exist([dir_fig,'probtrck_',meas,'_01_coreg_',Ori{end},'.png'],'file')
job1=[];
job1.subjID = subjID;
job1.base_query=[dir1,'${subj}/rt1w.nii'];
job1.over_query=[dir1,'${subj}/r',meas,'.nii'];
job1.seed_query=[dir1,'${subj}/rseed.nii'];
if ~exist([dir_fig,'probtrck_',meas,'_01_coreg_',Ori{end},'.png'],'file')
for k=1:3
job1.cfg=struct('slicedim',k, 'thres',0);
job1.fname_fig=[dir_fig,'probtrck_',meas,'_01_coreg_',Ori{k},'.png'];
imagethumbnail (job1);
end
end
end
end

%% fnirt
dir2=[JOB.dir_base,'/FA/'];
M=numel(IMG);
cd(dir1);
for n=1:N
subjid=subjID{n};
for j=1:M
fname_warp=[dir1,subjid,'/FA_to_target_warp.nii.gz'];
unix(['ln -sf ',dir2,'/',subjid,'_FA_to_target_warp.nii.gz ',fname_warp]);
fname_input =[dir1,subjid,'/r',IMG{j},'.nii'];
fname_output=[dir1,subjid,'/r',IMG{j},'_fnirted.nii.gz'];
if ~exist(fname_output,'file')
myunix(['$FSLDIR/bin/applywarp -i ',fname_input,' -o ',fname_output,...
' -r $FSLDIR/data/standard/FMRIB58_FA_1mm -w ',fname_warp]);
% now get a largest connected component
myfsl_lcc(fname_output, 0.1, fname_output);
end
end
end

% visual inspection!
if ~exist([dir_fig,'probtrck_',meas,'_02_fnirt_',Ori{end},'.png'],'file')
job1=[];
job1.subjID = subjID;
job1.base_query=[dir1,'/${subj}/rt1w_fnirted.nii.gz'];
job1.over_query=[dir1,'/${subj}/r',meas,'_fnirted.nii.gz'];
job1.seed_query=[dir1,'/${subj}/rseed_fnirted.nii.gz'];
for k=1:3
job1.cfg=struct('slicedim',k, 'thres',0);
job1.fname_fig=[dir_fig,'probtrck_',meas,'_02_fnirt_',Ori{k},'.png'];
imagethumbnail (job1);
end
end

% merge and mask (just whole brain)
dir_curr=pwd;
cd(dir3)
unix(['ln -sf ',dir_base,'/stats/mean_FA_mask.nii.gz ',dir3]);
unix(['fslmerge -t all_',meas,' ',dir1,'/*/r',meas,'_fnirted*']);
unix(['fslmaths all_',meas,' -mas mean_FA_mask all_',meas]);

%% skeletonise requires:

M=numel(MEAS);
for j=1:M
meas = MEAS{j};
unix(['ln -sf ',dir_base,'/stats/mean_FA.nii.gz ',dir3]);
unix(['ln -sf ',dir_base,'/stats/mean_FA_skeleton_mask_dst.nii.gz ',dir3]);
unix(['ln -sf ',dir_base,'/stats/all_FA.nii.gz ',dir3]);

if ~exist(['all_',meas,'_skeletonised.nii.gz'],'file')
unix(['tbss_skeleton -i mean_FA -p 0.2 mean_FA_skeleton_mask_dst ',...
'${FSLDIR}/data/standard/LowerCingulum_1mm ',...
'all_FA all_',meas,'_skeletonised -a all_',meas]);
unix(['fslmaths all_',meas,'_skeletonised -Tmean -bin ',...
'mean_',meas,'_skeletonised ']);

unix(['cp all_',meas,'.nii.gz ',dir_base,'/stats/']);
unix(['cp all_',meas,'_skeletonised.nii.gz ',dir_base,'/stats/']);
unix(['cp mean_',meas,'_skeletonised.nii.gz ',dir_base,'/stats/']);
end

% visual inspection!
if ~exist([dir3,'/',meas,'/',subjID{end},'.nii.gz'],'file')
[~,~]=mkdir([dir3,'/',meas,'_skeleton/']);
unix(['fslsplit all_',meas,'_skeletonised ./',meas,'_skeleton/ -t']);
for i=1:17
subjid=subjID{i};
movefile([meas,'_skeleton/',pad(i-1,4),'.nii.gz'], ...
[meas,'_skeleton/',subjid,'.nii.gz']);
end
end

if ~exist([dir_fig,'probtrck_',meas,'_03_skeleton_',Ori{end},'.png'],'file')
job1=[];
job1.subjID = subjID;
job1.base_query = [dir1,'/${subj}/rt1w_fnirted.nii.gz'];
job1.over_query = [dir3,'/',meas,'_skeleton/${subj}.nii.gz'];
job1.seed_query = [dir1,'/${subj}/rseed_fnirted.nii.gz'];
for k=1:3
job1.cfg=struct('slicedim',k, 'thres',0);
job1.fname_fig=[dir_fig,'probtrck_',meas,'_03_skeleton_',Ori{k},'.png'];
imagethumbnail (job1);
end
end
end

%% now GLM
numMeas = numel(MEAS);
rerun=ones(1,numMeas);
while ~~sum(rerun)
for k=1:numMeas
meas = MEAS{k};
[JOB.model_desc,JOB.COI] = fsss_model_desc(JOB.model,JOB.cidx);
% get links from "stats" directory to "GLM" directories
dir_glm = [dir_base,'/GLM/',meas,'_',JOB.glm_prefix,JOB.model_desc,'/'];
for s=1:2
X = char(JOB.model);
fname = [dir_glm,JOB.COI,'_tfce_corrp_tstat',num2str(s),'.nii.gz'];
if exist(fname,'file')
rerun(k) = 0;
end
end
end
for k=1:numMeas
meas = MEAS{k};
JOB.model_desc = fsss_model_desc(JOB.model,JOB.cidx);
% get links from "stats" directory to "GLM" directories
dir_glm = [dir_base,'/GLM/',meas,'_',JOB.glm_prefix,JOB.model_desc,'/'];
[~,~] = mkdir(dir_glm);
cd(dir_glm);
fnames=dir([dir_base,'/stats/*',meas,'*']); % all_FA.nii.gz, mean_FA.nii.gz, ...
for j=1:numel(fnames)
if ~fnames(j).isdir
unix(['ln -sf ',dir_base,'/stats/',fnames(j).name,' ',dir_glm,'.']);
end
end
JOB.fname_input = [dir_glm,'/all_',meas,'_skeletonised.nii.gz'];
ls(JOB.fname_input);

if isfield(JOB,'mask_thrs')
fn = [dir_glm,'/mean_',meas,'_skeletonised_mask',num2str(JOB.mask_thrs),'.nii.gz'];
unix(['fslmaths ',dir_glm,'/mean_',meas,'_skeletonised_mask.nii.gz ', ...
' -thr ',num2str(JOB.mask_thrs),' ',fn]);
JOB.fname_mask = fn;
else
JOB.fname_mask  = [dir_glm,'/mean_',meas,'_skeletonised_mask.nii.gz'];
end
ls(JOB.fname_mask);

JOB.dir_glm = dir_glm;
if rerun(k)
JOB.meas=meas;
JOB = myfsl_glm(JOB);
end
end
if ~~sum(rerun)
myunix(['waitForCONDORJobs.sh 60 ',JOB.qid]);
end
end

cd(dir_curr);
end

