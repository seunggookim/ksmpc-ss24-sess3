function Job = myants_combineTransforms(Job)
%
% antsApplyTransforms can be used to combine all Transforms and create a
% displacement field.
%
% eg: antsApplyTransforms -d 3 -o A2C.nii.gz -t A2B.h5 -t B2C.h5 -r C.nii
%
% this function just helps you to run the above line in MATLAB:
%
%
% Job.FnameOut = 'Warp_A2C.nii.gz'
% Job.Transforms = {'A2B.h5',useInverse; 'B2C.h5',useInverse};
% Job.FnameFixed = 'C.nii'
%
% antsApplyTransforms -d 3 -o [collapedWarp.nii.gz,1]
% -t [A2B.h5,0]
% -t [B2C.h5,0] -r mni_brain.nii.gz

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


%% check inputs
for itrans = 1:size(Job.Transforms,1)
  assert(isfile(Job.Transforms{itrans,1}), ...
    'Job.Transforms{%i,1}="%s" NOT FOUND!', itrans, ...
    Job.Transforms{itrans,1})
end
assert(isfile(Job.FnameFixed),...
  'Job.FnameFixed="%s" NOT FOUND!',Job.FnameFixed)

%% write command:
cmd = sprintf('antsApplyTransforms -d 3 -o [%s,1] -r %s ', ...
  Job.FnameOut, Job.FnameFixed);
for itrans = 1:size(Job.Transforms,1)
  cmd = [cmd sprintf(' -t [%s,%i] ', ...
    Job.Transforms{itrans,1}, Job.Transforms{itrans,2})];
end

%% run command:
system(cmd);
assert(isfile(Job.FnameOut))
end

function TEST()
cfg = struct('FnameOut',fn_reg_epi_to_mni, ...
  'FnameFixed','mni_funcref.nii.gz',...
  'Transforms',{{fn_reg_epi_to_t1w,0; fn_reg_t1w_to_mni,0}});

end
