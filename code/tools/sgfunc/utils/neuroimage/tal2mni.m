function xyz = tal2mni(xyz_tal)
% xyz_mni = tal2mni(xyz_tal)
%
% transformation matrix for pooled MNI-to-Talairach from:
% Lancaster et al. (2007). "Bias Between MNI and Talairach Coordinates
% Analyzed Using the ICBM-152 Brain Template". Human Brain Mapping
%
% (cc) 2015, code by sgKIM.

% Talairach to MNI
T=[
.9357 .0029 -.0072 -.0423
-.0065 .9396 -.0726 -1.3940
.0103 .0752 .8967 3.6475
0 0 0 1
];

if numel(xyz_tal) == 3
if size(xyz_tal,1) > size(xyz_tal,2)
xyz_tal=xyz_tal'; % 1 x 3
end
else
if size(xyz_tal,1) > size(xyz_tal,2)
xyz_tal=xyz_tal';
end
end

xyz = (T)*[(xyz_tal) ones(size(xyz_tal,1),1) ]';
xyz(4,:)=[];
xyz=xyz';
end
