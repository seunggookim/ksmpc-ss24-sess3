function JOB = myspm_rsfc(JOB)


JOB.subjID = fsss_subjID(JOB.subjID);

%% 1. fmri/t1w preprocessing: STC + unwarp/realignment (rigid-motion param)
%  .fwhm_mm
%  .fname_epi
%  .fname_t1w
%  .fname_vdm
JOB = myspm_fmriprep12_func(JOB);

%% 2. Scrubbing regressor by ART

JOB=[];
JOB.dir_base='/scr/vatikan3/APConn/rest12.410/';
JOB.subjID = subjID17;
JOB.name_epi='uarest410.nii';
JOB.name_rp ='rp_arest410.txt';
JOB.dir_figure='/scr/vatikan3/APConn/rest12.410/fig_denoising';
JOB.global_threshold=3;
JOB.motion_threshold=0.5;
JOB = myspm_art(JOB) % to make sure that you don't need scrubbing


%% 3. Compcor regressor
JOB.bpf1=[0 inf];
JOB=myy_compcor(JOB)


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
JOB = myspm_residual(JOB)

%% 5. registration [mni_2mm]



%% 6. smoothing [4 mm]



end
