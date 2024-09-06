function SL = mydsi_readSL(fname_dwi, fname_trk)
% SL = mydsi_readSL(fname_dwi, fname_trk)
%
% reads streamlines (SL) generated from DSI studio
% (dsi-studio.labsolver.org)
% 
% fname_dwi : the input NIFTI file
% fname_trk : streamlines exported in MATLAB format (.mat) with variables
%             "lengths" and "tracts"
%
% (cc) 2017, sgKIM.

load(fname_trk,'length','tracts')
nii = load_untouch_nii(fname_dwi); % is it always LAS? I guess so...?
d = size(nii.img);
T = [1 0 0 0; 0 -1 0 d(2)+1; 0 0 1 0; 0 0 0 1]; % LPS to LAS
tracts_xyz = ijk2xyz( ijk2uvw(tracts'+1, T), nii); % the voxel coordinates are zero-based

SL = [];
SL.num_sl = numel(length);
SL.num_pt_sl = length;
SL.length = zeros(1,SL.num_sl);
SL.coords = cell(1,SL.num_sl);
k = 1;
for j = 1:numel(SL.num_pt_sl)
 coords = [];
 for i = 1:SL.num_pt_sl(j)
  coords = [coords; tracts_xyz(k,:)];
  k = k+1;
 end
 SL.length(j) = sum(l2norm(diff(coords,1)));
 SL.coords{j} = coords;
end

end