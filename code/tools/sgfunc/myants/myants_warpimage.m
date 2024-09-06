function Job = myants_warpimage(Job)
% Job = myants_warpimage(Job)
%
% applies combined transformations using WarpImageMultiTransform.
%
% Forward warping means deformation from a moving image to a fixed image.
% Inverse warping means the opposite: from fixed to moving.
%
% Job (1x1) structure requires:
%  .FnameMoving
%  .FnameFixed
% (.fname_transforms) {'this_warp.nii.gz',isInv; 'that_warp.h5',isInv; ...}
% (.FnameOut)
% (.intp_order) = [0=NN, 1=3lin, 3=BSpline[3]]

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


%% check inputs:
Job.FnameMoving = getfullpath(Job.FnameMoving);
[p1,f1,e1] = myfileparts(Job.FnameMoving);
assert(isfile(Job.FnameMoving))

Job.FnameFixed = getfullpath(Job.FnameFixed);
[p2,f2,e2] = myfileparts(Job.FnameFixed);
assert(isfile(Job.FnameFixed))


if not(isfield(Job,'FnameOut'))
  Job.FnameOut = [p1,'/',f1,'_in_',f2,e1];
end

for itrans = 1:size(Job.transforms,1)
  assert(isfile(Job.transforms{itrans,1}))
end

cmd = sprintf('WarpImageMultiTransform 3 %s %s -R %s', ...
  Job.FnameMoving, Job.FnameOut, Job.FnameFixed);

if ~isfield(Job,'intp_order')
  Job.intp_order = 1;
end
switch Job.intp_order
  case 0
    cmd=[cmd ' --use-NN'];
  case 3
    cmd=[cmd ' --BSpline'];
  case 1
    % default without anyoption is trylinear
  otherwise
    error('BSpline[2] is not supported')
end
if ~isfield(Job,'supressfigure')
  Job.supressfigure = false;
end
if ~exist(Job.FnameOut,'file') || ~Job.supressfigure
  system(cmd);
  cfg = struct('contour',FnameFixed,...
    'fname_png',[p1,'/',f1,'_in_',f2,'.png']);
  slices(Job.FnameOut,[], cfg)
end
assert(isfile(Job.FnameOut))
% Job.FnameOut = FnameOut;
% Job.FnameWarp = FnameWarp;
% Job.FnameAff = FnameAff;

end
