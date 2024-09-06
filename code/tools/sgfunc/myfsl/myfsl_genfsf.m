fprintf('# FEAT version number\n');
fprintf('set fmri(version) 3.14\n');
fprintf('# Are we in MELODIC?\n');
fprintf('set fmri(inmelodic) 1\n');
fprintf('# Analysis level\n');
fprintf('# 1 : First-level analysis\n');
fprintf('# 2 : Higher-level analysis\n');
fprintf('set fmri(level) 1\n');
fprintf('# Which stages to run\n');
fprintf('# 0 : No first-level analysis (registration and/or group stats only)\n');
fprintf('# 7 : Full first-level analysis\n');
fprintf('# 1 : Pre-processing\n');
fprintf('# 2 : Statistics\n');
fprintf('set fmri(analysis) 7\n');
fprintf('# Use relative filenames\n');
fprintf('set fmri(relative_yn) 0\n');
fprintf('# Balloon help\n');
fprintf('set fmri(help_yn) 1\n');
fprintf('# Run Featwatcher\n');
fprintf('set fmri(featwatcher_yn) 1\n');
fprintf('# Cleanup first-level standard-space images\n');
fprintf('set fmri(sscleanup_yn) 0\n');
fprintf('# Output directory\n');
fprintf('set fmri(outputdir) "/scr/vatikan3/APConn/rest/1001"\n');
fprintf('# TR(s)\n');
fprintf('set fmri(tr) 1.5\n');
fprintf('# Total volumes\n');
fprintf('set fmri(npts) 420\n');
fprintf('# Delete volumes\n');
fprintf('set fmri(ndelete) 0\n');
fprintf('# Perfusion tag/control order\n');
fprintf('set fmri(tagfirst) 1\n');
fprintf('# Number of first-level analyses\n');
fprintf('set fmri(multiple) 1\n');
fprintf('# Higher-level input type\n');
fprintf('# 1 : Inputs are lower-level FEAT directories\n');
fprintf('# 2 : Inputs are cope images from FEAT directories\n');
fprintf('set fmri(inputtype) 2\n');
fprintf('# Carry out pre-stats processing?\n');
fprintf('set fmri(filtering_yn) 1\n');
fprintf('# Brain/background threshold, fprintf('set fmri(brain_thresh) 10\n');
fprintf('# Critical z for design efficiency calculation\n');
fprintf('set fmri(critical_z) 5.3\n');
fprintf('# Noise level\n');
fprintf('set fmri(noise) 0.66\n');
fprintf('# Noise AR(1)\n');
fprintf('set fmri(noisear) 0.34\n');
fprintf('# Motion correction\n');
fprintf('# 0 : None\n');
fprintf('# 1 : MCFLIRT\n');
fprintf('set fmri(mc) 1\n');
fprintf('# Spin-history (currently obsolete)\n');
fprintf('set fmri(sh_yn) 0\n');
fprintf('# B0 fieldmap unwarping?\n');
fprintf('set fmri(regunwarp_yn) 0\n');
fprintf('# EPI dwell time (ms)\n');
fprintf('set fmri(dwell) 0.7\n');
fprintf('# EPI TE (ms)\n');
fprintf('set fmri(te) 35\n');
fprintf('# fprintf('set fmri(signallossthresh) 10\n');
fprintf('# Unwarp direction\n');
fprintf('set fmri(unwarp_dir) y-\n');
fprintf('# Slice timing correction\n');
fprintf('# 0 : None\n');
fprintf('# 1 : Regular up (0, 1, 2, 3, ...)\n');
fprintf('# 2 : Regular down\n');
fprintf('# 3 : Use slice order file\n');
fprintf('# 4 : Use slice timings file\n');
fprintf('# 5 : Interleaved (0, 2, 4 ... 1, 3, 5 ... )\n');
fprintf('set fmri(st) 0\n');
fprintf('# Slice timings file\n');
fprintf('set fmri(st_file) ""\n');
fprintf('# BET brain extraction\n');
fprintf('set fmri(bet_yn) 1\n');
fprintf('# Spatial smoothing FWHM (mm)\n');
fprintf('set fmri(smooth) 5\n');
fprintf('# Intensity normalization\n');
fprintf('set fmri(norm_yn) 0\n');
fprintf('# Perfusion subtraction\n');
fprintf('set fmri(perfsub_yn) 0\n');
fprintf('# Highpass temporal filtering\n');
fprintf('set fmri(temphp_yn) 1\n');
fprintf('# Lowpass temporal filtering\n');
fprintf('set fmri(templp_yn) 0\n');
fprintf('# MELODIC ICA data exploration\n');
fprintf('set fmri(melodic_yn) 0\n');
fprintf('# Carry out main stats?\n');
fprintf('set fmri(stats_yn) 1\n');
fprintf('# Carry out prewhitening?\n');
fprintf('set fmri(prewhiten_yn) 1\n');
fprintf('# Add motion parameters to model\n');
fprintf('# 0 : No\n');
fprintf('# 1 : Yes\n');
fprintf('set fmri(motionevs) 0\n');
fprintf('set fmri(motionevsbeta) ""\n');
fprintf('set fmri(scriptevsbeta) ""\n');
fprintf('# Robust outlier detection in FLAME?\n');
fprintf('set fmri(robust_yn) 0\n');
fprintf('# Higher-level modelling\n');
fprintf('# 3 : Fixed effects\n');
fprintf('# 0 : Mixed Effects: Simple OLS\n');
fprintf('# 2 : Mixed Effects: FLAME 1\n');
fprintf('# 1 : Mixed Effects: FLAME 1+2\n');
fprintf('set fmri(mixed_yn) 2\n');
fprintf('# Number of EVs\n');
fprintf('set fmri(evs_orig) 1\n');
fprintf('set fmri(evs_real) 1\n');
fprintf('set fmri(evs_vox) 0\n');
fprintf('# Number of contrasts\n');
fprintf('set fmri(ncon_orig) 1\n');
fprintf('set fmri(ncon_real) 1\n');
fprintf('# Number of F-tests\n');
fprintf('set fmri(nftests_orig) 0\n');
fprintf('set fmri(nftests_real) 0\n');
fprintf('# Add constant column to design matrix? (obsolete)\n');
fprintf('set fmri(constcol) 0\n');
fprintf('# Carry out post-stats steps?\n');
fprintf('set fmri(poststats_yn) 1\n');
fprintf('# Pre-threshold masking?\n');
fprintf('set fmri(threshmask) ""\n');
fprintf('# Thresholding\n');
fprintf('# 0 : None\n');
fprintf('# 1 : Uncorrected\n');
fprintf('# 2 : Voxel\n');
fprintf('# 3 : Cluster\n');
fprintf('set fmri(thresh) 3\n');
fprintf('# P threshold\n');
fprintf('set fmri(prob_thresh) 0.05\n');
fprintf('# Z threshold\n');
fprintf('set fmri(z_thresh) 2.3\n');
fprintf('# Z min/max for colour rendering\n');
fprintf('# 0 : Use actual Z min/max\n');
fprintf('# 1 : Use preset Z min/max\n');
fprintf('set fmri(zdisplay) 0\n');
fprintf('# Z min in colour rendering\n');
fprintf('set fmri(zmin) 2\n');
fprintf('# Z max in colour rendering\n');
fprintf('set fmri(zmax) 8\n');
fprintf('# Colour rendering type\n');
fprintf('# 0 : Solid blobs\n');
fprintf('# 1 : Transparent blobs\n');
fprintf('set fmri(rendertype) 1\n');
fprintf('# Background image for higher-level stats overlays\n');
fprintf('# 1 : Mean highres\n');
fprintf('# 2 : First highres\n');
fprintf('# 3 : Mean functional\n');
fprintf('# 4 : First functional\n');
fprintf('# 5 : Standard space template\n');
fprintf('set fmri(bgimage) 1\n');
fprintf('# Create time series plots\n');
fprintf('set fmri(tsplot_yn) 1\n');
fprintf('# Registration to initial structural\n');
fprintf('set fmri(reginitial_highres_yn) 0\n');
fprintf('# Search space for registration to initial structural\n');
fprintf('# 0   : No search\n');
fprintf('# 90  : Normal search\n');
fprintf('# 180 : Full search\n');
fprintf('set fmri(reginitial_highres_search) 90\n');
fprintf('# Degrees of Freedom for registration to initial structural\n');
fprintf('set fmri(reginitial_highres_dof) 3\n');
fprintf('# Registration to main structural\n');
fprintf('set fmri(reghighres_yn) 1\n');
fprintf('# Search space for registration to main structural\n');
fprintf('# 0   : No search\n');
fprintf('# 90  : Normal search\n');
fprintf('# 180 : Full search\n');
fprintf('set fmri(reghighres_search) 90\n');
fprintf('# Degrees of Freedom for registration to main structural\n');
fprintf('set fmri(reghighres_dof) BBR\n');
fprintf('# Registration to standard image?\n');
fprintf('set fmri(regstandard_yn) 1\n');
fprintf('# Use alternate reference images?\n');
fprintf('set fmri(alternateReference_yn) 0\n');
fprintf('# Standard image\n');
fprintf('set fmri(regstandard) "/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain"\n');
fprintf('# Search space for registration to standard space\n');
fprintf('# 0   : No search\n');
fprintf('# 90  : Normal search\n');
fprintf('# 180 : Full search\n');
fprintf('set fmri(regstandard_search) 90\n');
fprintf('# Degrees of Freedom for registration to standard space\n');
fprintf('set fmri(regstandard_dof) 12\n');
fprintf('# Do nonlinear registration from structural to standard space?\n');
fprintf('set fmri(regstandard_nonlinear_yn) 0\n');
fprintf('# Control nonlinear warp field resolution\n');
fprintf('set fmri(regstandard_nonlinear_warpres) 10 \n');
fprintf('# High pass filter cutoff\n');
fprintf('set fmri(paradigm_hp) 288\n');
fprintf('# Total voxels\n');
fprintf('set fmri(totalVoxels) 208158720\n');
fprintf('# Number of lower-level copes feeding into higher-level analysis\n');
fprintf('set fmri(ncopeinputs) 0\n');
fprintf('# 4D AVW data or FEAT directory (1)\n');
fprintf('set feat_files(1) "/scr/vatikan3/APConn/rest/1001/aurest"\n');
fprintf('# Add confound EVs text file\n');
fprintf('set fmri(confoundevs) 0\n');
fprintf('# Subject's structural image for analysis 1\n');
fprintf('set highres_files(1) "/scr/vatikan3/APConn/rest/1001/t1w_brain"\n');
fprintf('# Resampling resolution\n');
fprintf('set fmri(regstandard_res) 4\n');
fprintf('# Variance-normalise timecourses\n');
fprintf('set fmri(varnorm) 1\n');
fprintf('# Automatic dimensionality estimation\n');
fprintf('set fmri(dim_yn) 1\n');
fprintf('# Output components\n');
fprintf('set fmri(dim) 1\n');
fprintf('# 1 : Single-session ICA\n');
fprintf('# 2 : Multi-session temporal concatenation\n');
fprintf('# 3 : Multi-session tensor TICA\n');
fprintf('set fmri(icaopt) 1\n');
fprintf('# Threshold IC maps\n');
fprintf('set fmri(thresh_yn) 1\n');
fprintf('# Mixture model threshold\n');
fprintf('set fmri(mmthresh) 0.5\n');
fprintf('# Output full stats folder\n');
fprintf('set fmri(ostats) 0\n');
fprintf('# Timeseries and subject models\n');
fprintf('set fmri(ts_model_mat) ""\n');
fprintf('set fmri(ts_model_con) ""\n');
fprintf('set fmri(subject_model_mat) ""\n');
fprintf('set fmri(subject_model_con) ""\n');
fprintf('##########################################################\n');
fprintf('# Now options that don't appear in the GUI\n');
fprintf('# Alternative (to BETting) mask image\n');
fprintf('set fmri(alternative_mask) ""\n');
fprintf('# Initial structural space registration initialisation transform\n');
fprintf('set fmri(init_initial_highres) ""\n');
fprintf('# Structural space registration initialisation transform\n');
fprintf('set fmri(init_highres) ""\n');
fprintf('# Standard space registration initialisation transform\n');
fprintf('set fmri(init_standard) ""\n');
fprintf('# For full FEAT analysis: overwrite existing .feat output dir?\n');
fprintf('set fmri(overwrite_yn) 0\n');