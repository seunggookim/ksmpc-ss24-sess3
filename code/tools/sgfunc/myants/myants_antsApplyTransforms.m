function Job = myants_antsApplyTransforms(Job)
% job = myants_antsApplyTransforms(job)
%
% job (1x1) structure requires:
%  .FnameFixed
%  .FnameMoving
% (.FnameOut)
% (.odt)
% (.interpolation)   'linear' | 'NearestNeighbor' | {'BSpline[<order=3>]'}
%                    | 'LanczosWindowedSinc' | and more...
%  .transforms       {Nx2} filenames in 1st col, useInverse in 2nd col
% (.IsCreateFig)     0 = no figure | 1 = warped slices (default)
%
% NOTE: for transform A -> B -> C: {A-to-B},{B-to-C}
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
assert(isfile(Job.FnameFixed))
assert(isfile(Job.FnameMoving))

%% Set up defaults
if ~isfield(Job,'FnameOut')
  [p1,f1,e1] = myfileparts(Job.FnameFixed);
  [p2,f2,e2] = myfileparts(Job.FnameMoving);
  prefix = [p2,'/',f2,'_in_',f1];
  Job.FnameOut = [prefix,'.nii.gz'];
end
if ~isfield(Job,'interpolation')
  Job.interpolation = 'BSpline[3]';
end
if ~isfield(Job,'IsCreateFig')
  Job.IsCreateFig = true;
end

%% Data type:
if ~isfield(Job,'odt')
  info = niftiinfo(Job.FnameMoving);
  switch (info.Datatype)
    case 'int16'
      datatype = 'short';
    case {'single', 'int32'}
      datatype = 'float';
    otherwise
      error('FIXME!')
  end
else
  datatype = Job.odt;
end

%% antsApplyTransforms
cmd = sprintf('antsApplyTransforms -d %i -i %s -r %s -o %s -n %s -u %s ', ...
  3, Job.FnameMoving, Job.FnameFixed, Job.FnameOut, Job.interpolation, datatype);
% concatenate transforms:
for i=1:size(Job.Transforms,1)
  cmd = [cmd, sprintf(' --transform [%s,%i]', Job.Transforms{i,1}, Job.Transforms{i,2})];
end
if not(Job.IsCreateFig)
  cmd = [cmd, ' -v 0 '];
end
system(cmd);

%% visualize results
if Job.IsCreateFig
  [p3,f3,~] = fileparts(Job.FnameOut);
  cfg = struct('fname_png',[p3,'/',f3,'.png'],'contour',Job.FnameFixed);
  slices(Job.FnameOut,[], cfg)
end

end
