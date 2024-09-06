function JOB = myfsl_probtrackx2 (JOB)
% JOB = myfsl_probtrackx2 (JOB)
%
% JOB requires:
%  .dname_bedpostx
%  .name_probtrackx
% (.dname_probtrackx)
%  .fname_seed
%  .fname_target
% (.fname_waypoint)
% (.fname_stop)
% (.optarg) any optional arguments for probtrackx2
%
% Waypoint ROI: If a curve does not go through, it is discarded.
% Exclusion ROI: If a curve goes through, it is discarded.
% Termination ROI: If a curve goes through, it is terminated.
%
% NETWORK MODES:
%
% OUTPUT MATRIX
% (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FDT/UserGuide#BEDPOSTX)
% Matrix1: SEED-to-SEED conn matrix
%          The seed mask(s) can represent grey matter,
%          so this would be all GM to GM connectivity.
% Matrix2: SEED-to-TARGET2
%          The seed can be a grey matter region, and mask2 the rest of the
%          brain (can be low res). The results can then be used for blind
%          classification of the seed mask.
% Matrix3: TARGET3-to-TARGET3 or TARGET3-to-LRTARGET3
%          The 'target' masks can be the whole of grey matter, and the seed
%          mask white matter. This option allows more sensitive
%          reconstruction of grey-to-grey connectivity as pathways are
%          seeded from all their constituting locations, rather than just
%          their end-points as in Matrix1.
% Matrix4: TARGET4-to-SEED or TARGET4-to-COLMASK4
%          DtiMask-to-Seed (Oxford Sparse Format?)
%
%
% (cc) 2015-2019, sgKIM. solleo@gmail.com  https://ggooo.wordpress.com

%% check data:
[p1,f1,e1] = myfileparts(JOB.dname_dwi);
dn_bedpostx = [p1,'/',f1,e1,'.bedpostX'];
ls([dn_bedpostx,'/merged_f1samples.nii.gz']);
% dn_probtrackx = [p1,'/',seedname];
if ~isfield(JOB,'dname_probtrackx')
  dn_probtrackx = [p1,'/',f1,e1,'.prbtrk.',JOB.name_probtrackx];
else
  dn_probtrackx = JOB.dname_probtrackx;
end
[~,~] = mkdir(dn_probtrackx);

if ~iscell(JOB.fname_seed)
  [p1,f1,e1] = myfileparts(JOB.fname_seed);
  fn_seed = [p1,'/',f1,e1];
  ls(fn_seed)
else
  fid = fopen([dn_probtrackx,'/seedslist.txt'], 'w');
  for k=1:numel(JOB.fname_seed)
    [p1,f1,e1] = myfileparts(JOB.fname_seed{k});
    fn_seed = [p1,'/',f1,e1];
    ls(fn_seed)
    fprintf(fid,'%s\n', fn_seed);
  end
  fclose(fid);
  fn_seed = [dn_probtrackx,'/seedslist.txt'];
end
fn_mask = [dn_bedpostx,'/nodif_brain_mask.nii.gz'];


%%
setenv('FSLOUTPUTTYPE','NIFTI_GZ')
cmd=fullfile(getenv('FSLDIR'),'bin','probtrackx2');
% inputs:
cmd=[cmd,' --samples=',dn_bedpostx,'/merged',...
  ' --mask=',fn_mask,' --seed=',fn_seed];

% outputs:
cmd=[cmd,' --dir=',dn_probtrackx,' --forcedir',...; % output directory name
  ' --opd'];%,... % output path distribution
%   ' --pd']; % Correct path distribution for the length of the pathways:
% the connectivity distribution is the expected length of the pathways that
% cross each voxel times the number of samples that cross it


cmd=[cmd,' --ompl'];%,... % Output mean path length from seed 
% writes the mean of the lengths of all the accepted streamlines in mm. 
% It will generate extra files (with suffix _lengths) for all the 
% distribution maps and all the connectivity matrices generated.

%   ' --os2t']; % Output seeds to targets

% 	' --fopd',... % Other mask for binning tract distribution

% connectivity matrices:


% streamline generation & thresholding control:
if ~isfield(JOB,'nsamples')
  JOB.nsamples = 5000;
end
cmd=[cmd,' --nsamples=',num2str(JOB.nsamples),... % # of streamlines per voxel {5000}
  ' --nsteps=2000',... % # of steps per sample
  ' --steplength=0.5',... % step length in mm
  ' --distthresh=0.0',... % distance threshold in mm
  ' --cthr=0.2']; % curvature threshold
%    ' --usef']; % use anisotropy to constrain tracking

% waypoint mask: ONLY include a streamline when it hits this mask
if isfield(JOB,'fname_waypoint')
  % create 'waypoint mask list' file
  fid = fopen([dn_probtrackx,'/waypointslist.txt'], 'w');
  if ~iscell(JOB.fname_waypoint)
    JOB.fname_waypoint = {JOB.fname_waypoint};
  end
  for k=1:numel(JOB.fname_waypoint)
    [p1,f1,e1] = myfileparts(JOB.fname_waypoint{k});
    fn_wpt = [p1,'/',f1,e1];
    fprintf(fid,'%s\n', fn_wpt);
  end
  fclose(fid);
  cmd= [cmd,' --waypoints=',dn_probtrackx,'/waypointslist.txt',...
    ' --waycond=AND'];
end

% exclusion mask: exclude a streamline when it hits this mask



% stop mask: terminate tracking when it hits this mask
if isfield(JOB,'fname_stop')
  [p1,f1,e1] = myfileparts(JOB.fname_stop);
  fn_stp = [p1,'/',f1,e1];
  cmd = [cmd,' --stop=',fn_stp];
end

cmd=[cmd,' --verbose=1 '];
if isfield(JOB,'optarg')
  cmd=[cmd,' ',JOB.optarg];
end
if ~exist([dn_probtrackx,'/fdt_paths.nii.gz'],'file')
unix(cmd);
end

pwd0 = pwd;
cd(dn_probtrackx)
pathdist = load_untouch_nii('fdt_paths.nii.gz');
brain = load_untouch_nii([JOB.dname_dwi,'/bmnodiff_unwarped.nii']);
lengths = load_untouch_nii('fdt_paths_lengths.nii.gz');

cfg=struct('contournum',3, 'contourcolor',[.5 .5 .5],...
  'colormap',hot, 'caxisprc',[0 100]);

cfg.colorbartitle='# tracks [smp]';
cfg.fname_png = 'numTracks.png';
slices(pathdist.img,cfg,brain.img);

cfg.colorbartitle='mean dist [mm]';
cfg.colormap=[0 0 0; parula];
cfg.fname_png = 'meanDist.png';
slices(lengths.img,cfg,brain.img);

cfg.colorbartitle='# tracks x mean dist [smp.mm]';
cfg.colormap=hot;
cfg.fname_png = 'numTracksByDist.png';
slices(pathdist.img.*lengths.img,cfg,brain.img);

if ~isempty(dir('particle*'))
  mkdir manyparticles
  unix('mv particle* manyparticles/')
end
cd(pwd0)

JOB.dname_probtrackx = dn_probtrackx;
JOB.fname_numtracks = [dn_probtrackx,'/fdt_paths.nii.gz'];
JOB.fname_lengths = [dn_probtrackx,'/fdt_paths_lengths.nii.gz'];

end
