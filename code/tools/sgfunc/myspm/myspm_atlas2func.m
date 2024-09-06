%function JOB = myspm_atlas2func (JOB)

subjid='2001';
JOB.name_t1w='Brain.nii';
JOB.name_epi='uarest410.nii';
JOB.dir_base='/scr/vatikan3/APConn/rest12.410/';
JOB.prefix_ants='/scr/vatikan3/APConn/ants_t1w_n17/t1w';
JOB.fname_atlas='/scr/vatikan1/skim/matlab/conn/rois/atlas.nii';

path1=[fullfile(JOB.dir_base,subjid),'/'];
JOB.fname_t1w=[path1,JOB.name_t1w];
[~,name1,~]=fileparts(JOB.fname_t1w);
JOB.fname_epi=[path1,'mean',JOB.name_epi];
[~,name3,~]=fileparts(JOB.fname_atlas);

%%

% 1. 
% atlas@mni => @t1w
job1=[];
job1.fname_ref = JOB.fname_t1w;
job1.fname_in  = JOB.fname_atlas;
job1.fname_out = [path1,name3,'_in_',name1,'.nii'];
job1.xfm(1).prefix = [JOB.prefix_ants,subjid];
job1.xfm(1).isInverse = 1; % MNI to t1w
job1.intp='NearestNeighbor';
myants_applyxfm(job1);

% atlas@t1w => @func
job1=[];
job1.name_fixed  = JOB.fname_epi;
job1.name_moving = JOB.fname_t1w;
job1.name_others = [path1,name3,'_in_',name1,'.nii'];
job1.interp=0;
myspm_coreg(job1);

% mask it with GM>0.5
JOB.fname_atlas = [path1,'o',name3,'_in_',name1,'.nii'];
gm = load_uns_nii([path1,'oc1t1w.nii']); % in integer format!
gm.img = zeroone(gm.img);
atlas = load_uns_nii(JOB.fname_atlas);
atlas.img = uint16((gm.img>0.5) .* atlas.img);
atlas.hdr.dime.datatype=512;
save_untouch_nii(atlas, JOB.fname_atlas);

%%
% %% But not a good idea...? There are some errors, but also in MNI space too.
% % it should be better to warp it, and mask with GM>0.50
% job1.name_fixed  = JOB.fname_t1w;
% job1.name_moving = JOB.fname_epi;
% job1.prefix='g';
% myspm_coreg(job1)
% job1=[];
% job1.name_def='y_t1w.nii';
% job1.name_moving='gmeanuarest410.nii';
% job1.vox_mm=[2 2 2];
% tic
% myspm_deform1(job1)
% 
