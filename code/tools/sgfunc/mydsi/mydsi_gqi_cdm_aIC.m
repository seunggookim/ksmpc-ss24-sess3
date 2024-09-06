function job = mydsi_gqi_cdm_aIC(job)
% job = mydsi_gqi_cdm_aIC(job)
%
% job requires:
%  .dir_subj
% (.name_data)
% (.mdsr) diffusion sampling length ratio. For Zhang's data 0.5 used
% (.name_seed) default='aIC' from /Volumes/macExJ_4TB/VV/DSI/data/roi
% (.name_roi) default='VV_frntstr11' from /Volumes/macExJ_4TB/VV/DSI/data/roi
% (.num_tracks_k) # of tracks in k (x1000) (default=100)
% (.tip) topology-informed pruning iteration (default=0)

dir_subj = job.dir_subj;
if ~isfield(job,'mdsr'), job.mdsr='0.5'; end
mdsr = job.mdsr;
if isnumeric(mdsr), mdsr=num2str(mdsr); end
if ~isfield(job,'name_seed'), job.name_seed = 'aIC'; end
name_seed = job.name_seed;
if ~isfield(job,'dn_roi')
  dn_roi = [getenv('HOME'),'/Dropbox/DSI/data/roi'];
end
if ~isfield(job,'name_roi'), job.name_roi = 'VV_frntstr11'; end
name_roi = job.name_roi;
if ~isfield(job,'num_tracks_k'), job.num_tracks_k = 100; end
if ~isfield(job,'tip'), job.tip=0; end
if ~isfield(job,'name_data'), job.name_data='dsi171_loca'; end

dir0=pwd;
cd (dir_subj)
% dsicmd='srun -n2 --mem 16g --time 01:00:00 --pty --x11 singularity run /applications/singularity/images/dsi-studio-latest.img ';
if ismac
  dsicmd='/Applications/dsi_studio.app/Contents/MacOS/dsi_studio ';
elseif isunix
  dsicmd='/usr/local/dsi_studio_64/dsi_studio ';
end
% dsicmd='/Applications/dsi_studio.app/Contents/MacOS/dsi_studio ';


%% Check inputs:
fn_data = [job.name_data,'.nii.gz'];
if exist([job.name_data,'.nii'],'file') ...
    && exist([job.name_data,'.nii.gz'],'file')
  unix(['gzip ',job.name_data,'.nii'])
end
ls(fn_data)
if ~exist('bmnodiff.nii','file') && exist('bmnodiff.nii.gz','file')
  unix('gunzip bmnodiff.nii.gz')
end
ls('bmnodiff.nii')
ls('y_nodiff.nii')


%% Check seed & roi files:
fn_seed=['i',name_seed,'.nii'];% aIC.nii';
fn_roi=['i',name_roi,'.nii'];%VV_frntstr11.nii';
if ~exist(fn_seed,'file')
  fn_seed_mni = [dn_roi,'/',name_seed,'.nii'];
  myspm_denorm(struct('fname_fixed','bmnodiff.nii',...
    'fname_deform','y_nodiff.nii',...
    'fname_moving',fn_seed_mni,'interp',0,'savepwd',1));
end
if ~exist(fn_roi,'file')
  fn_roi_mni = [dn_roi,'/',name_roi,'.nii'];
  myspm_denorm(struct('fname_fixed','bmnodiff.nii',...
    'fname_deform','y_nodiff.nii',...
    'fname_moving',fn_roi_mni,'interp',0,'savepwd',1));
  fn_roi_txt = [dn_roi,'/',name_roi,'.txt'];
  copyfile(fn_roi_txt,[pwd,'/'])
end


%% Import data
if ~exist([job.name_data,'.src.gz'],'file')
  unix([dsicmd,' --action=src ',...
    '--source=',fn_data,' --bval=',job.name_data,'.bval ',...
    '--bvec=',job.name_data,'.bvec ',...
    '--output=',job.name_data,'.src.gz'])
end
setenv('FSLOUTPUTTYPE','NIFTI_GZ')
unix(['fslmaths bmnodiff.nii.gz -bin -dilF -ero nodiff_mask.nii.gz ',...
  '-odt char'])


%% Recon: method=4 (GQI)
[~,fn_fib]=myls(...
  [job.name_data,'.src.gz.odf8.f5.*fy.rdi.gqi.',mdsr,'.fib.gz']);
if isempty(fn_fib)
  unix([dsicmd,' --action=rec ',...
    '--source=',job.name_data,'.src.gz --mask=nodiff_mask.nii.gz ',...
    '--method=4 ',... % 4=GQI 7=QSDR (Q-Space diffeomorphic reconstruction (QSDR) [11] is the generalization of GQI that allows users to construct spin distribution functions (SDFs, a kind of ODF) in any given template space (e.g. MNI space). )
    '--param0=',mdsr,' ',... % For QSDR, param0 is the mean diffusion distance ratio, and param1 is the output resolution in mm.
    '--delete_repeat=1 ',...
    '--odf_order=8 --num_fiber=5'])
end


%% Tractography: wholebrain for sanity check
%  Euler with topology-informed pruning
%  (Runge-Kutta4 was more sensitive to the crazy artifacts...)
%  auto-gen method: "A diffusion spectrum imaging scheme was used, and a total of 256 diffusion sampling were acquired. The maximum b-value was 0 s/mm2. The in-plane resolution was 2.4 mm. The slice thickness was 2.4 mm. The diffusion data were reconstructed in the MNI space using q-space diffeomorphic reconstruction (Yeh et al., Neuroimage, 58(1):91-9, 2011) to obtain the spin distribution function (Yeh et al., IEEE TMI, ;29(9):1626-35, 2010).  A diffusion sampling length ratio of 1.25 was used The restricted diffusion was quantified using restricted diffusion imaging (Yeh et al., MRM, 77:603–612 (2017)). A deterministic fiber tracking algorithm (Yeh et al., PLoS ONE 8(11): e80713) was used. A seeding region was placed at whole brain. The angular threshold was 30 degrees. The step size was 0.05 mm. The anisotropy threshold was 0.01. The fiber trajectories were smoothed by averaging the propagation direction with a percentage of the previous direction. The percentage was randomly selected from 0% to 95%. Tracks with length shorter than 30 or longer than 300 mm were discarded. A total of 10000 tracts were calculated. Topology-informed pruning (Yeh et al. 2018) was applied to the tractography with 3 iteration(s) to remove false connections. parameter_id=0AD7233C9A99193FD7B35D3FCDCC4C3Db803FbF041b96431027cb01cb03d"

[~,fn_fib]=myls([job.name_data,'.src.gz.odf8.f5*rdi.gqi.',mdsr,'.fib.gz']);
if isempty(fn_fib)
  error('recon result not found in %s',pwd)
end
fn_fib=fn_fib{1};
fn_trk=[fn_fib,'.100k.wholebrain.trk.gz'];
if ~exist([fn_trk,'.2.png'],'file')
  unix([dsicmd,' --action=trk ',...
    '--source=',fn_fib,' ',...
    '--fiber_count=100000 --fa_threshold=0.02 ',...
    '--turning_angle=30 --step_size=.5 --smoothing=0.5 --min_length=30 ',...
    '--max_length=300 --tracking_method=0 --tip_iteration=3 ',...
    '--output=',fn_trk])
  unix([dsicmd,' --action=vis ',...
    '--source=',fn_fib,' --track=',fn_trk,' ',...
    '--cmd="',...
    'set_view,0;save_image,',fn_trk,'.0.png;',...
    'set_view,1;set_view,1;save_image,',fn_trk,'.1.png;',...
    'set_view,2;save_image,',fn_trk,'.2.png"']);
end


%% Tractography: anterior IC
% fn_trk=[fn_fib,'.',nk.aIC.trk.gz'];
fn_trk=sprintf('%s.%ik.tip%i.%s.trk.gz', fn_fib, job.num_tracks_k, job.tip, job.name_seed);
if ~exist([fn_trk,'.',job.name_roi,'.count.end.connectivity.mat'],'file')
  unix([dsicmd,' --action=trk ',...
    ' --source=',fn_fib,' --seed=',fn_seed,...
    ' --fiber_count=',num2str(job.num_tracks_k*1000),...
    ' --fa_threshold=0.02 ',...
    ' --turning_angle=50 --step_size=.5 --smoothing=0.5',...
    ' --min_length=30 ',...
    ' --max_length=300 --tracking_method=0',...
    ' --tip_iteration=',num2str(job.tip),...
    ' --output=',fn_trk,...
    ' --connectivity=',fn_roi,' ',.... %VV_frntstr11 ',...
    ' --connectivity_value=qa,count,ncount --connectivity_type=end'])
  unix([dsicmd,' --action=vis ',...
    '--source=',fn_fib,' --add_surface=',fn_seed,' --track=',fn_trk,' ',...
    '--cmd="',...
    'set_view,0;save_image,',fn_trk,'.0.png;',...
    'set_view,1;set_view,1;save_image,',fn_trk,'.1.png;',...
    'set_view,2;save_image,',fn_trk,'.2.png"']);
end


%%
cd(dir0)
end
