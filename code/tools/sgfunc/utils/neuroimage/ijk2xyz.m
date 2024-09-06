function Xyz = ijk2xyz(Ijk1, Info)
% converts 1-based voxel-subscripts to world-coordinates (mm)
%
% SYNTAX
% Xyz = ijk2xyz(Ijk1, Info)
%
% INPUTS
% Ijk1   [3xN]  MATLAB's 1-based voxel-subscripts of N points
% Info   (1x1)  a structure from NIFTIINFO
%      | '1xN'  a filename of the NIFTI file
%
% OUTPUT
% Xyz    [3xN]  World-coordinates (mm) of N points
%
% (cc) 2022, dr.seunggoo.kim@gmail.com

if ischar(Info)
  Info = niftiinfo(Info);
else
  error('Info can be either a structure or a string')
end
Ijk0 = Ijk1 - 1;
Xyz = [Ijk0 ones(size(Ijk0,1),1)]*Info.Transform.T;
Xyz(:,4) = [];

end


%%
function TEST()
FnMni = fullfile(getenv('FSLDIR'),'data','standard','MNI152_T1_1mm.nii.gz');
xyz = ijk2xyz([90 126 72]+1, FnMni);
assert(isequal(xyz, [0 0 0])) % world origin
xyz = ijk2xyz([0 0 0]+1, FnMni); % voxel origin
assert(isequal(xyz, [90 -126 -72])) 
end
