function y = myspm_readvox (fname, xyz_mm, ijk_flag)
%y = myspm_readvox (fname, xyz_mm, ijk_flag)
% (cc) 2016, sgKIM.
if exist('ijk_flag')&&ijk_flag
%ijk = reshape(xyz_mm,[3 1]);
ijk = xyz_mm';
else
ijk = xyz2ijk(xyz_mm, fname)';
end
V = spm_vol(fname);
y=spm_get_data(V,round(ijk));
end

% %% Sanity check:
% function amIinsane
% clear
% n=1
% load /scr/vatikan4/conmus3/mat/info23.mat
% dir1=['/scr/vatikan4/conmus3/GLM/',mrtID{n}];
% fname=[dir1,'/swudata.nii'];
% 
% nii=load_uns_nii(fname);
% xyz=[12 -49.5 -5];
% ijk=xyz2ijk(xyz,nii) % 27 26 19
% y = myspm_readvox(fname, xyz);
% 
% y1=squeeze(nii.img(27,26,19,:));
% sum(y-double(y1)) % zero. no mistake in myspm_readvox
% end
