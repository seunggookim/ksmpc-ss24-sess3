function FD_j = myy_FD_jenkinson(fname_rp, fname_mni)
% FD_j = myy_FD_jenkinson(fname_rp, fname_mni)
%
% computes frame-wise displacement according to the Jenkinson.
%
% Reference: Jenkinson, M., Bannister, P., Brady, M., Smith, S., 2002. Improved optimization for the robust and accurate linear registration and motion correction of brain images. Neuroimage 17, 825-841.
% Note: The lines below are from y_FD_Jenkinson.m by Yan Chao-Gen
% (cc?) 2016, sgKIM.

rmax = 80.0; %The default radius (as in FSL) of a sphere represents the brain

RP=load(fname_rp);
nTimePoint=size(RP,1);
sinq1=sin(RP(:,4));
sinq2=sin(RP(:,5));
sinq3=sin(RP(:,6));
cosq1=cos(RP(:,4));
cosq2=cos(RP(:,5));
cosq3=cos(RP(:,6));

if ~exist('fname_mni','var')
fname_mni=[getenv('FSLDIR'),'/data/standard/MNI152_T1_1mm.nii.gz'];
hdr = load_untouch_header_only(fname_mni);
center = [hdr.hist.srow_x; hdr.hist.srow_y; hdr.hist.srow_z; 1 1 1 1]*[(hdr.dime.dim(2:4)')./2; 1];
center = center(1:3);
end

FD_J = zeros(nTimePoint-1,1);
for t=2:nTimePoint
M_RigidBodyTransform_0 = RigidBodyTransform(cosq1,sinq1,cosq2,sinq2,cosq3,sinq3,t,RP);
M_RigidBodyTransform_1 = RigidBodyTransform(cosq1,sinq1,cosq2,sinq2,cosq3,sinq3,t,RP);  

MA1=(M_RigidBodyTransform_1);
MA2=(M_RigidBodyTransform_0);

M = MA1*inv(MA2) - eye(4);

A = M(1:3,1:3);

T = M(1:3,4);

FD_J(t-1) = sqrt(rmax*rmax/5*trace(A'*A) + (T+A*center)'*(T+A*center));

end
FD_J=[0;FD_J]; %The FD_Jenkinson at time point t means the movement from time point t-1 to time point t. (Put the FD_Jenkinson for the first time point to "0".)

end

function M_rigid = RigidBodyTransform(cosq1,sinq1,cosq2,sinq2,cosq3,sinq3,t,RP)
M1=[1       0        0     0;...
0    cosq1(t)  sinq1(t)  0;...
0    -sinq1(t) cosq1(t)  0;...
0       0        0     1;];

M2=[cosq2(t)  0    sinq2(t)     0;...
0        1       0        0;...
-sinq2(t) 0    cosq2(t)     0;...
0       0        0        1;];

M3=[cosq3(t)   sinq3(t)   0     0;...
-sinq3(t)  cosq3(t)   0     0;...
0           0       1     0;...
0           0       0     1;];

MT=[1    0     0     RP(t,1);...
0    1     0     RP(t,2);...
0    0     1     RP(t,3);...
0    0     0     1;];

M_RigidBodyTransform=MT*M1*M2*M3;
end
