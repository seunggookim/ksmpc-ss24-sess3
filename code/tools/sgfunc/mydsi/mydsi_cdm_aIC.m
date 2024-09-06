function job = mydsi_cdm_aIC(job)
% dir_ = job.dir_sub;
if ismac
  dsicmd = '/Applications/dsi_studio.app/Contents/MacOS/dsi_studio';
end
fn_dsi = job.fname_dsi;
[dir_data,prefix_dsi,~] = myfileparts(fn_dsi);
dir0 = pwd;
cd (dir_data)


%% Import data
if ~exist('dsi.src.gz','file')
  unix([dsicmd,' --action=src --source=',fn_dsi,...
    ' --bval=',prefix_dsi,'.bval --bvec=',prefix_dsi,'.bvec',...
    ' --output=',prefix_dsi,'.src.gz'])
end


%% Create no diffusion mask
unix(['fslroi ',prefix_dsi,'.nii nodiff.nii 0 1'])
if ~exist('bmnodiff.nii','file')
  myspm_seg12(struct('fname_t1w','nodiff.nii'),'ss')
end
unix(['fslmaths bmnodiff.nii -bin -dilF -ero nodiff_mask.nii -odt char'])


%% Recon: CDM
[~,fn_fib]=myls('dsi.src.gz.odf8.f5*cdm.qsdr.0.6.R*.fib.gz');
if isempty(fn_fib)
  unix(['/Applications/dsi_studio.app/Contents/MacOS/dsi_studio --action=rec ',...
    '--source=dsi.src.gz --mask=nodiff_mask.nii --method=7 ',...
    '--param0=0.6 --reg_method=3 --delete_repeat=1',...
    '--odf_order=8 --num_fiber=5 --half_sphere=1'])
end
%% Tractography: Euler with topology-informed pruning
%  (Runge-Kutta4 was more sensitive to the crazy artifacts...)
%{
 A diffusion spectrum imaging scheme was used, and a total of 256 diffusion sampling were acquired. The maximum b-value was 0 s/mm2. The in-plane resolution was 2.4 mm. The slice thickness was 2.4 mm. The diffusion data were reconstructed in the MNI space using q-space diffeomorphic reconstruction (Yeh et al., Neuroimage, 58(1):91-9, 2011) to obtain the spin distribution function (Yeh et al., IEEE TMI, ;29(9):1626-35, 2010).  A diffusion sampling length ratio of 1.25 was used The restricted diffusion was quantified using restricted diffusion imaging (Yeh et al., MRM, 77:603–612 (2017)). A deterministic fiber tracking algorithm (Yeh et al., PLoS ONE 8(11): e80713) was used. A seeding region was placed at whole brain. The angular threshold was 30 degrees. The step size was 0.05 mm. The anisotropy threshold was 0.01. The fiber trajectories were smoothed by averaging the propagation direction with a percentage of the previous direction. The percentage was randomly selected from 0% to 95%. Tracks with length shorter than 30 or longer than 300 mm were discarded. A total of 10000 tracts were calculated. Topology-informed pruning (Yeh et al. 2018) was applied to the tractography with 3 iteration(s) to remove false connections. parameter_id=0AD7233C9A99193FD7B35D3FCDCC4C3Db803FbF041b96431027cb01cb03d
%}
[~,fn_fib]=myls('dsi.src.gz.odf8.f5*cdm.qsdr.0.6.R*.fib.gz');
fn_fib=fn_fib{1};
fn_seed='/Volumes/macExJ_4TB/VV/DSI/roi/ALIC.nii';
fn_trk='cdm.0.6.100k.IC.trk.gz';
if ~exist(fn_trk,'file')
  unix(['/Applications/dsi_studio.app/Contents/MacOS/dsi_studio --action=trk ',...
    '--source=',fn_fib,' --seed=',fn_seed,' ',...
    '--fiber_count=100000 --export=fa,statistics --fa_threshold=0.02 ',...
    '--turning_angle=50 --step_size=.5 --smoothing=0.5 --min_length=30 ',...
    '--max_length=300 --tracking_method=0 --tip_iteration=3 ',...
    '--output=',fn_trk])
end
%'--cluster=2,12,2.4,EMclus.txt ',...

% unix(['/Applications/dsi_studio.app/Contents/MacOS/dsi_studio --action=trk ',...
% '--parameter_id=0AD7233C9A99193FD7B35D3FCDCC4C3Db803FbF041b96431027cb01cb03d  ',...
% '--source=',fn_src,' ',...
% '--seed=JHU-WhiteMatter-labels-1mm:Anterior_limb_of_internal_capsule_R ',...
% '--output=cdm.IC_R.trk.gz'])

%% visualize:
fn_seed='/Volumes/macExJ_4TB/VV/DSI/roi/ALIC.nii';
if ~exist([fn_trk,'.2.png'],'file')
unix(['/Applications/dsi_studio.app/Contents/MacOS/dsi_studio --action=vis ',...
  '--source=',fn_fib,' ',...
  '--add_surface=',fn_seed,' ',...
  '--track=',fn_trk,' ',...
  '--cmd="',...
  'set_view,0;save_image,',fn_trk,'.0.png;',...
  'set_view,1;set_view,1;save_image,',fn_trk,'.1.png;',...
  'set_view,2;save_image,',fn_trk,'.2.png"']);
end
%% connectivity matrix:
unix(['/Applications/dsi_studio.app/Contents/MacOS/dsi_studio --action=trk ',...
'--source=',fn_fib,' --tract=',fn_trk,' --connectivity=VV_frntstr11 ',...
'--connectivity_value=qa,count,ncount --connectivity_type=end'])

%%
cd(dir0)
end
