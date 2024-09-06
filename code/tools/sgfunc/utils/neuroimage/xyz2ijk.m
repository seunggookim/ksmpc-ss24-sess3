function Ijk1 = xyz2ijk(Xyz, Info)
% converts world-coordinates (mm) to 1-based voxel-subscripts
% 
% SYNTAX
% Ijk1 = xyz2ijk(Xyz, Info)
%
% INPUTS
% Xyz    [3xN]  World-coordinates (mm) of N points
% Info   (1x1)  a structure from NIFTIINFO
%      | '1xN'  a filename of the NIFTI file
%
% OUTPUT
% Ijk1   [3xN]  MATLAB's 1-based voxel-subscripts of N points
%
% (cc) 2022, dr.seunggoo.kim@gmail.com

if ischar(Info)
  Info = niftiinfo(Info);
else
  error('Info can be either a structure or a string')
end
Ijk0 = [Xyz ones(size(Xyz,1),1)]/Info.Transform.T;
Ijk1 = Ijk0(:,1:3) + 1;

end


%%
function TEST()
FnMni = fullfile(getenv('FSLDIR'),'data','standard','MNI152_T1_1mm.nii.gz');
ijk1 = xyz2ijk([0 0 0], FnMni); % world origin
assert(isequal(ijk1, [90 126 72]+1))
ijk1 = xyz2ijk([90 -126 -72], FnMni)
assert(isequal(ijk1, [0 0 0]+1)); % voxel origin
end
