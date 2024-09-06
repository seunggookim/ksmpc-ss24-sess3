function JOB = myspm_rsfc_example(JOB)


JOB.subjID = fsss_subjID(JOB.subjID);

%% 1. fmri/t1w preprocessing: STC + unwarp/realignment (rigid-motion param)
%  .fwhm_mm
%  .fname_epi
%  .fname_t1w
%  .fname_vdm
JOB = myspm_fmriprep12_func(JOB); % 'a' for STC; 'u' for unwarping/realign

%% 2. Scrubbing regressor by ART

JOB=[];
JOB.dir_base='/scr/vatikan3/APConn/rest12.410/';
JOB.subjID = subjID17;
JOB.name_epi='uarest410.nii';
JOB.name_rp ='rp_arest410.txt';
JOB.dir_figure='/scr/vatikan3/APConn/rest12.410/fig_denoising';
JOB.global_threshold=3;
JOB.motion_threshold=0.5;
JOB = myspm_art(JOB); % to make sure that you don't need scrubbing

%% 3. Compcor regressor
JOB.bpf1=[0 inf];
JOB=myy_compcor(JOB);

%% 4. compare results from various regressions

JOB.param_cc='wmcsf99_n16d1v1b0.00-Inf';
JOB.name_cc=['cc_',JOB.param_cc,'_eigenvec.txt'];
JOB.param_art='3.0std_0.5mm';
JOB.name_art=['art_regression_outliers_and_movement_uarest410_',JOB.param_art,'.mat'];
JOB.name_rp='rp_arest410.txt';
JOB.fname_gmmask='oc1t1w_99.nii';
JOB.bpf2=[0.01 0.10];
JOB.covset=[1 2 3 4];
JOB.cov_idx=4; % i'll go with global-signal!
JOB = myspm_denoise(JOB);  % 'r' for residual; 'f' for filtering

%% 5. bring parcellation into function

JOB=[];
JOB.dir_base = '/scr/vatikan3/APConn/rest12.410/';
JOB.dir_fs   = '/scr/vatikan3/APConn/FSspm12';
JOB.name_t1w = 'Brain.nii';
JOB.name_epi = 'meanuarest410.nii';
JOB.subjID = subjID;
myspm_aparc(JOB)

%% 6. compute correlation without smoothing
% smoothing (doesn't seem to be a good idea for ROIs. Averaging already
% cancels out tiny noise. Even 2-voxel-fhwm increases skewness and drops
% gof from ks-test.... (but i'll use it for seed-based correlation!@
JOB=[];
JOB.dir_base = '/scr/vatikan3/APConn/rest12.410/';
JOB.subjID = subjID;
JOB.dir_figure = '/scr/vatikan3/APConn/figures.local/crrmtx_sg';
JOB.name_res='fruarest410.nii';
myspm_crrmtx(JOB)




%% also FS-projection...


%% 7. register only 3D image! 
% (registering 4-D double images TAKES TOO TOO LONG! 10 mins per subject)
% coreg 2.3mm,410scans,single func => strc (2mm): takes 3 min
% deform




end
