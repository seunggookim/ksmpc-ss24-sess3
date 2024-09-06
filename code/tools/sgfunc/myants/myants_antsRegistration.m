function JOB = myants_antsRegistration(JOB)
% JOB = myants_antsRegistration(JOB)
%
% JOB requires:
%  .fname_fixed
%  .fname_moving
% (.dname_out)        [default: path of fname_moving]
% (.reg_stages)       0=rigid, 1=+affine, 2=+SyN [default]
% (.interpolation)    'linear' | 'NearestNeighbor' | 'BSpline[<order=3>]' | 
%                     'LanczosWindowedSinc' | and more...
% (.reg_SyNtransform) 'SyN[0.1,3,0]' for inter-subject reg [default]
%                     'SyN[0.1,3,0.5]' for intra-subject reg
%                     
%                     
% REF: https://github.com/ANTsX/ANTs/wiki/Anatomy-of-an-antsRegistration-call
%
% (cc) 2019, sgKIM

%% check LD_LIBRARY_PATH
%{
  NOTE: Starting up MATLAB will put MATLAB's lib paths before ANTs lib path.
  If MATLAB has its own ITK libraries (e.g., HPC) this can mask ANTs's lib.
  So I just need to switch the order in that case:
%}
LD_LIBRARY_PATH = getenv('LD_LIBRARY_PATH');
PATH_ = strsplit(LD_LIBRARY_PATH,':');
ind = contains(upper(PATH_),'ANTS');
assert(sum(ind), 'No ANTs library in LD_LIBRARY_PATH?!');
PATH_ = [PATH_(ind), PATH_(~ind)];
LD_LIBRARY_PATH = cell2mat(strcat(PATH_,':'));
setenv('LD_LIBRARY_PATH', LD_LIBRARY_PATH(1:end-1));


%% Check inputs
[p1,f1,e1] = myfileparts(JOB.fname_fixed);
[p2,f2,e2] = myfileparts(JOB.fname_moving);
fn_fixed=[p1,'/',f1,e1];
assert(isfile(fn_fixed))
fn_moving=[p2,'/',f2,e2];
assert(isfile(fn_moving))
% hdr_fixed = load_untouch_header_only(fn_fixed);
Info = niftiinfogz(fn_fixed);
% voxsize_fixed = prod(hdr_fixed.dime.pixdim(2:4));
voxsize_fixed = prod(Info.PixelDimensions(1:3));
% hdr_moving = load_untouch_header_only(fn_moving);
Info = niftiinfogz(fn_moving);
% voxsize_moving = prod(hdr_moving.dime.pixdim(2:4));
voxsize_moving = prod(Info.PixelDimensions(1:3));
if voxsize_fixed > voxsize_moving
  warning(['fixedImage is in a lower resolution than movingImage. ',...
    'Because shrink factor is based on fixedImage, this may cause a coarse ',...
    'registration. LowRes-to-highRes registration is recommended.']);
end
if ~isfield(JOB,'dname_out')
  JOB.dname_out = p2;
end
if ~isfield(JOB,'reg_stages')
  JOB.reg_stages = 2; % 0=rigid, 1=+affine, 2=+SyN [default]
end
prefix = [JOB.dname_out,'/',...
  f2,'_to_',f1,'_stage',num2str(JOB.reg_stages),'_'];
fn_reg=[prefix,'Composite.h5'];

fn_fwdwarped = [prefix,'fwdwarped.nii.gz'];
fn_invwarped = [prefix,'invwarped.nii.gz'];


%% antsRegistration
cmd=['antsRegistration ']; % --dimensionality 3 
if ~isfield(JOB,'dimensionality')
  JOB.dimensionality = 3;
end
cmd=[cmd,' --dimensionality ',num2str(JOB.dimensionality)];
% output options
cmd=[cmd,...
  ' --write-composite-transform 1',... % this is composite (0th to 3rd)
  ' --float 0 ',...
  ' --output [',prefix,',',fn_fwdwarped,',',fn_invwarped,'] '];
if ~isfield(JOB,'interpolation')
  JOB.interpolation = 'linear';
end
cmd=[cmd,' --interpolation ',JOB.interpolation,' '];

% preprocessing (clipping, intensity matching, coordinate initialization)
if ~isfield(JOB,'useHistogramMatching')
  JOB.useHistogramMatching = 0; % 0=for inter-modal, 1=for intra-modal reg
end
cmd=[cmd,...
  ' --winsorize-image-intensities [0.005,0.995] ',... % clipping outliers (%)
  ' --use-histogram-matching ',num2str(JOB.useHistogramMatching),' ']; 
cmd=[cmd, ' --initial-moving-transform [',fn_fixed,',',fn_moving,',1] '];
% matching 0=midpoint of bounding boxes, 1=center of mass, 2=origin

% 0th transform: rigid
if JOB.reg_stages >= 0
  cmd=[cmd,...
    ' --transform Rigid[0.1] ',... % gradient (0.1-0.25 recommended)
    ' --metric MI[',fn_fixed,',',fn_moving,',1,64,Regular,0.25] '];
  % metric[fixed, moving, weight, bins, sampling, samplingPercentage]
  %   weight=1 for inter-modal reg; 0.7 for T1-to-T1; 0.3 for T2-to-T2
  cmd=[cmd, ' --convergence [1000x500x250x100,1e-6,10] '];
  % [maxIterByLevel,threshold,convergenceWindowSize]
  % iteration stops if the metric improves less than threshold for the last
  %   convergenceWindowSize iterations or it reaches maxIterByLevel
  cmd=[cmd, ' --shrink-factors 8x4x2x1 '];
  % subsampling factors based on FIXED image, thus always register LOW-to-HIGH,
  %   not high-to-low!
  cmd=[cmd, ' --smoothing-sigmas 3x2x1x0vox '];
  % c.f. fwhm = 2.35*sigma
end

% 1st transform: affine
if JOB.reg_stages >= 1
  cmd=[cmd,... % parameters are same as rigid
    ' --transform Affine[0.1] ',...
    ' --metric MI[' fn_fixed ',' fn_moving ',1,32,Regular,0.25] ',...
    ' --convergence [1000x500x250x100,1e-6,10] ',...
    ' --shrink-factors 8x4x2x1 ',...
    ' --smoothing-sigmas 3x2x1x0vox '];
end

% 2nd transform: symmetric image normalization (SyN)
if JOB.reg_stages >= 2
  if ~isfield(JOB,'reg_SyNtransform')
    JOB.reg_SyNtransform = 'SyN[0.1,3,0]';
  end
  cmd=[cmd,' --transform ',JOB.reg_SyNtransform,' '];
  % [gradientStep,updateFieldVarianceInVoxelSpace,totalFieldVarianceInVoxelSpace]
  % gradient step - allowed movement at each iteration (0.1-0.25
  % recommended) updateFieldVarianceInVoxelSpace -  a smoothing penalty on
  % updated (this iteration) gradient field by neighbouring voxels (# of
  % neighbours) totalFieldVarianceInVoxelSpace - a smoothing penalty on
  % total (from the beginning to the current iteration) gradient field
  cmd=[cmd,' --metric CC[' fn_fixed ',' fn_moving ',1,4] '];
  % CC[fixedImage,movingImage,weight,radius]
  % cross-correlation between spheres around each voxel in two images
  cmd=[cmd,... % parameters are similar to rigid/affine
    ' --convergence [100x70x50x20,1e-6,10] ',...
    ' --shrink-factors 8x4x2x1 ',...
    ' --smoothing-sigmas 3x2x1x0vox'];
end
cmd = [cmd,' --verbose 2>&1 | tee ',prefix,'.log'];
if not(isfile(fn_reg))
  tic
  fid = fopen([prefix,'cmd.sh'],'w');
  fprintf(fid,'%s\n',cmd);
  fclose(fid);
  fprintf('[%s] START (%s)\n', mfilename, datestr(now,31));
  disp(cmd);
  system(cmd);
  fprintf('[%s] END (%s) ', mfilename, datestr(now,31));
  toc
end
JOB.fname_reg = fn_reg;
JOB.fname_fwdwarped = fn_fwdwarped;
JOB.fname_invwarped = fn_invwarped;

%% visualize results
[p3,f3,~] = myfileparts(fn_fwdwarped);
cfg = struct('fname_png',[p3,'/',f3,'.png'],'contour',fn_fixed);
slices(fn_fwdwarped,[], cfg)
[p3,f3,~] = myfileparts(fn_invwarped);
cfg = struct('fname_png',[p3,'/',f3,'.png'],'contour',fn_moving);
slices(fn_invwarped,[], cfg)
end
