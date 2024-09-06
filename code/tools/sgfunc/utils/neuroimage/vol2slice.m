function [Vi,u,w, IJK] = vol2slice(V, T, ras, method, forceras)
% [Vi,u,w] = vol2slice(V, T, ras, [method], forceras)
%
% [INPUT]
% V is a 3-D volume in I x J x K
%
% RAS is a 3-D RAS(right/anterior/superior) coordinate but only one
% dimension should have non-nan value (e.g, [nan nan 0]) for 1-based voxel
% coordinates.
%
% T is a 4x4 transformatrix such that [x,y,z,1]' = T*[i,j,k,1]'
% where i,j,k are 1-based indices in volume V
%
% METHOD is a string for interpolation option for INTERPN: 'nearest'
% 'linear' 'spline' 'cubic' 'makima'
% or maximum intensity projection 'mip'
%
% FORCERAS is a logical flag to make the Vi in RAS coordinate or not
%
% [OUTPUT]
% Vi is actually W-by-U not (U-by-W) because IMAGESC shows your
% row-by-column matrix on the [-Y,+X] space (ie. axis ij). But on the
% [+X,+Y] space (axis xy), a sane image should be in column-by-row (C
% convension). It kills you when you think about it. So just forget all
% about them until next time you have to handle this. After all, you just
% want to draw some pretty pictures (on a "sane" cooridnate system).
%
% For now, you can just use this output like:
%     imagesc(u.axis, w.axis, Vi); axis xy;
%     xlabel(u.axisname); ylabel(w.axisname);
% for a nice slice of a brain volume, which works well with other MATLAB
% functions to add annotations!
%
% U and W have .axis (coordinates in mm) and .axisname
%
% (cc) 2020, sgKIM. solleo@gmail.com

% 2021-05-14: BUG:one voxel misaligned.. =+=
%


R = inv(T); % World (RAS?LPI?...) to Voxel matrix
[rasdim2ijkdim] = findmaxrowpercol(abs(R(1:3,1:3)));
[ijkdim2rasdim] = findmaxrowpercol(abs(T(1:3,1:3))); % ijk2ras

rasd = ras;
rasd(isnan(ras)) = 0;
ijk = R*[rasd ones(size(ras,1),1)]'; %#ok<MINV> % ras2ijk
rasdim = find(~isnan(ras));
ijkdim = rasdim2ijkdim(rasdim); %#ok<FNDSB>

% RAS-coordinates to IJK-axes in voxel-space
IJK = cell(1,3); Xi = []; axisnames = 'RAS';
for i = 1:3
  IJK{i} = 1:size(V,i);
  G = ones(size(V,i),4);
  G(:,i) = IJK{i};
  U = T*G'; % ijk2ras
  Xi(i).axis = U(ijkdim2rasdim(i),:);
  Xi(i).axisname = axisnames(ijkdim2rasdim(i));
end

if ~exist('method','var') && isempty(method)
  method = 'linear';
end

if ~exist('forceras','var')
  forceras = true;
end

% Make it a singleton for the dimension we want to squeeze:
IJK{ijkdim} = ijk(ijkdim);
Vi = getslice(V, IJK, method, ijkdim);
%{
Image is YX... not XY, so Vi is transposed. Now the column-major
convension in C makes sense (for XY images; but not for matrices). Now I
want to say IMAGESC is crazy, so you have to transpose 2-D images
everytime you use IMAGESC. But...
%}
Vi = squeeze(Vi)';

% other dimensions
leftdims = setdiff(1:3, ijkdim);
u = Xi(leftdims(1));
w = Xi(leftdims(2));

% Here we force RAS direction on the 2-D slices again:
if forceras
  if strcmp(u.axisname,'S') || ...% rotated sagittal or rotated coronal
      (strcmp(u.axisname,'A') && strcmp(u.axisname,'R')) % rotated axial
    Vi = Vi';
    u = Xi(leftdims(2));
    w = Xi(leftdims(1));
  end
end

end

function Vijk = getslice(V, IJK, method, ijkdim)
if strcmp(method,'mip') % maximum intensity projection
  Vijk = max(V, [], ijkdim);
else
  Vijk = interpn(double(V), IJK{1}, IJK{2}, IJK{3}, method);
end
end

function rows = findmaxrowpercol(X)
rows = nan(size(X,1),1);
for icol = 1:size(X,2)
  [~,rows(icol)] = max(X(:,icol));
end
end

function TEST()
%%
vox2ras = eye(4);
vox2ras(1:3,4) = -1; % now 1-based: (1,1,1)th -> (0,0,0)mm
vol = zeros([5 5 5]);
vol(3,3,3) = 10;  % (3,3,3)th -> (2,2,2)mm
[Vi,u,w] = vol2slice(vol, vox2ras, [nan 2 nan], 'nearest');
assert(Vi(3,3) == 10)

vox2ras = eye(4);
vox2ras(1:3,4) = -1; % now 1-based: (1,1,1)th -> (0,10,0)mm
vox2ras(2,4) = vox2ras(2,4) + 10;
vol = zeros([5 5 5]);
vol(3,3,3) = 10;   % now 1-based: (3,3,3)th -> (2,12,2)mm
[Vi,u,w] = vol2slice(vol, vox2ras, [nan 12 nan], 'nearest');
assert(Vi(3,3) == 10)

%% BUT IS THE WORLD-COORDINATE OF THE FIRST VOXEL REALLY [0,0,0]mm?

% in FSL: (90,-126,-72)
% in FS: (90,-126,-72)

info = niftiinfo('~/MNI152_T1_1mm.nii');
T = info.Transform.T';
vol = niftiread(info);
[Vi,u,w,ijk] = vol2slice(vol, vox2ras_0to1(T), [nan 0 nan], 'nearest')
% U: from +90 to -91
% W: from -72 to 109
imagesc(u.axis,w.axis,Vi); axis xy image;
assert(Vi(w.axis==0, u.axis==0) == 3912)

%%
t10p5 = '/usr/local/fsl/data/standard/MNI152_T1_0.5mm.nii.gz';
t12 = '/usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz';
slices(t10p5, t12, struct('facealpha',0.5,'colormap',hot))

%%
func = '/home/sgk/github/musencfmri/data/spm/sub-01/func/happy/s6wuaepi_happy_skip5.nii';
t1w = '/home/sgk/github/musencfmri/data/spm/sub-01/anat/wbmt1w.nii';
slices(t1w, func, struct('facealpha',0.7,'colormap',hot))

%%
func = '/home/sgk/github/musencfmri/data/spm/sub-01/func/happy/s6wuaepi_happy_skip5.nii';
t1w = '/home/sgk/Documents/MATLAB/spm12/canonical/avg152T1.nii';
slices(t1w, func, struct('facealpha',0.7,'colormap',hot))
%%
t1w = '/home/sgk/github/musencfmri/data/spm/sub-01/anat/wbmt1w.nii';
mni = '/home/sgk/Documents/MATLAB/spm12/canonical/avg152T1.nii';
slices(mni, t1w, struct('facealpha',0.7,'colormap',hot,'thres',100))

%%
subj = '/home/sgk/github/musencfmri/data/spmants/sub-01/bmt1w_to_mni_brain_stage2_fwdwarped.nii.gz';
mni = '/home/sgk/Documents/MATLAB/spm12/canonical/avg152T1.nii';
slices(mni, subj, struct('facealpha',0.7,'colormap',hot))
%%
resliced = '/home/sgk/github/musencfmri/data/spm/sub-01/func/happy/wuaepi_happy_skip5.nii';
mni = '/home/sgk/Documents/MATLAB/spm12/canonical/avg152T1.nii';
slices(mni, resliced, struct('facealpha',0.7,'colormap',hot))
%%
resliced = '/home/sgk/github/musencfmri/data/spm/sub-01/func/happy/wuaepi_happy_skip5.nii';
mni = '/home/sgk/Documents/MATLAB/spm12/canonical/avg152T1.nii';
slices(resliced, [], struct('facealpha',0.7,'colormap',hot, 'contour',mni))
%%
resliced = '/home/sgk/github/musencfmri/data/spm/sub-01/func/happy/wuaepi_happy_skip5.nii';
mni = '/home/sgk/Documents/MATLAB/spm12/canonical/avg152T1.nii';
slices(mni, [], struct('facealpha',0.7,'colormap',hot, 'contour',resliced))

%%

end


%{
%% EXAMPLES:














%%
mri = load_nifti('/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz')

thal = load_nifti('/usr/local/fsl/data/atlases/Thalamus/Thalamus-prob-1mm.nii.gz')'

clf
colormap(gray)
ax = axeslayout([1 6],'tight','tight');
for i = 1:6
  y = -35.5+6*i;
[Vi,u,w] = vol2slice(mri.vol, mri.vox2ras, [nan y nan],'lineer');
axespos(ax,i)
imagesc(u.axis,w.axis,Vi); axis xy image;
set(gca,'xtick',[],'ytick',[])
% xlabel(u.axisname); ylabel(w.axisname);
% axis([-20 20 20 60])
xlabel(sprintf('Y = %.0f',y)); grid off;

axespos(ax,i)
[Qi,u,w] = vol2slice(thal.vol(:,:,:,1), thal.vox2ras, [nan y nan],'linear');
Qi(~Qi) = nan;
h = pcolor(u.axis, w.axis, Qi);
h.LineStyle = 'none';
colormap(gca,turbo)
axis off image

end

%%
mri = load_nifti('/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz')

thal = load_nifti('/usr/local/fsl/data/atlases/Thalamus/Thalamus-maxprob-thr50-2mm.nii.gz')'

clf
colormap(gray)
ax = axeslayout([1 6],'tight','tight');
for i = 1:6
  y = -35.5+6*i;
[Vi,u,w] = vol2slice(mri.vol, mri.vox2ras, [nan y nan],'linaer');
axespos(ax,i)
imagesc(u.axis,w.axis,Vi); axis xy image;
set(gca,'xtick',[],'ytick',[])
% xlabel(u.axisname); ylabel(w.axisname);
% axis([-20 20 20 60])
xlabel(sprintf('Y = %.0f',y)); grid off;
axis([-40 40 -40 40])

axespos(ax,i)
[Qi,u,w] = vol2slice(thal.vol, thal.vox2ras, [nan y nan],'nearest');
Qi(~Qi) = nan;
h = pcolor(u.axis, w.axis, Qi);
h.LineStyle = 'none';
colormap(gca,brewermap(numel(unique(thal.vol(:))),'dark2'))
axis off image
axis([-40 40 -40 40])

end

%%

mri = load_nifti('/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz');
thal = load_nifti('/usr/local/fsl/data/atlases/Thalamus/Thalamus-maxprob-thr50-2mm.nii.gz');
prob = load_nifti('/usr/local/fsl/data/atlases/Thalamus/Thalamus-prob-1mm.nii.gz');

clf
colormap(gray)
ax = axeslayout([1 6],'tight','tight');
for i = 1:6
  y = -35.5+6*i;
[Vi,u,w] = vol2slice(mri.vol, mri.vox2ras, [nan y nan],'linear');
axespos(ax,i)
imagesc(u.axis,w.axis,Vi); axis xy image;
set(gca,'xtick',[],'ytick',[])
xlabel(sprintf('Y = %.0f',y)); grid off;
axis([-40 40 -40 40])


axespos(ax,i)
[Qi,u,w] = vol2slice(prob.vol(:,:,:,2), prob.vox2ras, [nan y nan],'linear');
Qi(~Qi) = nan;
h = pcolor(u.axis, w.axis, Qi);
h.LineStyle = 'none';
colormap(gca,flipud(brewermap(256,'spectral')))
axis off image
axis([-40 40 -40 40])

hold on
[Qi,u,w] = vol2slice(thal.vol, thal.vox2ras, [nan y nan],'nearest');
[~,h] = contour(u.axis, w.axis, convn(~~Qi,ones(3),'same'), 1);
% [~,h] = contour(u.axis, w.axis, Qi>0, 1);
h.LineWidth = 1;
h.LineColor = 'k';


end

%%
mri = load_nifti('/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz');
thal = load_nifti('/usr/local/fsl/data/atlases/Thalamus/Thalamus-maxprob-thr50-2mm.nii.gz');
prob = load_nifti('/usr/local/fsl/data/atlases/Thalamus/Thalamus-prob-1mm.nii.gz');

clf
colormap(gray)
ax = axeslayout([1 6],'tight','tight');
for i = 1:6
  y = -35.5+6*i;
[Vi,u,w] = vol2slice(mri.vol, mri.vox2ras, [nan y nan],'linear');
axespos(ax,i)
imagesc(u.axis,w.axis,Vi); axis xy image;
set(gca,'xtick',[],'ytick',[])
xlabel(sprintf('Y = %.0f',y)); grid off;
axis([-40 40 -40 40])

axespos(ax,i)
[Qi,u,w] = vol2slice(prob.vol(:,:,:,2), prob.vox2ras, [nan y nan],'linear');
Qi(~Qi) = nan;
h = pcolor(u.axis, w.axis, Qi);
h.LineStyle = 'none';
colormap(gca,flipud(brewermap(256,'spectral')))
axis off image
axis([-40 40 -40 40])

hold on
[Qi,u,w] = vol2slice(thal.vol, thal.vox2ras, [nan y nan],'nearest');
[~,h] = contour(u.axis, w.axis, convn(~~Qi,ones(3),'same'), 1);
% [~,h] = contour(u.axis, w.axis, Qi>0, 1);
h.LineWidth = 1;
h.LineColor = 'k';
end


%% "Glass brain"------------------------------------------------------
mri = load_nifti('~/MNI152_T1_1mm.nii');
thal = load_nifti('/usr/local/fsl/data/atlases/Thalamus/Thalamus-maxprob-thr50-2mm.nii.gz');
prob = load_nifti('/usr/local/fsl/data/atlases/Thalamus/Thalamus-prob-1mm.nii.gz');

clf
coords = [0 nan nan; nan -10 nan; nan nan 5];

ax = axeslayout([1 3],'tight','tight');
for i = 1:3
  
  axespos(ax,i)
  [Qi,u,w] = vol2slice(prob.vol(:,:,:,2), prob.vox2ras, coords(i,:),'mip');
  Qi(~Qi) = nan;
  h = pcolor(u.axis, w.axis, Qi);
  h.LineStyle = 'none';
  h.FaceAlpha = 0.9;
  colormap(gca, (brewermap(256,'Blues')))
  grid on
  set(gca,'xticklabel',[],'yticklabel',[])
  axis xy image
  
  axespos(ax,i)
  [Vi,u,w] = vol2slice(mri.vol, mri.vox2ras, coords(i,:),'linear');
  [~,h] = contour(u.axis, w.axis, convn(Vi,ones(3),'same'), 1);
  h.LineColor = [.5 .5 .5];
  set(gca,'color','none')
  axis xy image; box on;
  set(gca,'xtick',[],'ytick',[])
  xyzdim = find(~isnan(coords(i,:)));
  xyzname = 'XYZ';
  xlabel(sprintf('%s = %.0f',xyzname(xyzdim), coords(i,xyzdim)));
  grid on
  
  hold on
  [Qi,u,w] = vol2slice(thal.vol, thal.vox2ras, coords(i,:),'nearest');
  [~,h] = contour(u.axis, w.axis, convn(~~Qi,ones(3),'same'), 1);
  % [~,h] = contour(u.axis, w.axis, Qi>0, 1);
  h.LineWidth = 0.5;
  h.LineColor = [.8 .5 0];
  
  
end


%}