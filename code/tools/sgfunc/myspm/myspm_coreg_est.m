function out = myspm_coreg_est(fname_moving, fname_fixed, fname_dat)
% T = myspm_coreg_est(fname_moving, fname_fixed, fname_dat)
%
% returns a rigid-body transformation matrix from [fname_moving] to [fname_fixed]
% [fname_dat] can be optionally defined.
%
% (cc) 2017, sgKIM.

if ~nargin,  help myspm_coreg_est;  return; end
[p1,f1,~]=fileparts(fname_moving);
[~,f2,~]=fileparts(fname_fixed);
if ~exist('fname_dat','var')
 fname_dat=[p1,'/coregmat_',f1,'_to_',f2,'.dat'];
end
if ~exist(fname_dat,'file')
 VG=spm_vol(fname_fixed);
 VF=spm_vol(fname_moving);
 x = spm_coreg(VG,VF);
 T = inv(VF.mat\spm_matrix(x(:)')*VG.mat);
 save(fname_dat,'T','-ascii');
else
 T=load(fname_dat,'T','-ascii');
end
out=[];
out.T=T;
out.fname_dat = fname_dat;
end
