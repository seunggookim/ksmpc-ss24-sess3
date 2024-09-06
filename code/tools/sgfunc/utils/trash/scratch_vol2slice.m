

% Make a slice with WORLD-coordinates

T = mni.vox2ras1;
V = mni.vol;
x = [nan nan 0]; % position x in WORLD-coordinate


[I,J,K] = ndgrid(1:size(V,1), 1:size(V,2), 1:size(V,3));
XYZ = (T*[I(:),J(:),K(:),ones(numel(V),1)]')';
X = reshape(XYZ(:,1),size(V));
Y = reshape(XYZ(:,2),size(V));
Z = reshape(XYZ(:,3),size(V));

R = inv(T);
[ijkdim2rasdim,~] = find(R(:,1:3));



%%
% figure;
% h = slice(I,J,K,V,[1.5],[10.90],[30.5],'linear');
% h(1).LineStyle='none';
% h(2).LineStyle='none';
% h(3).LineStyle='none';
% xlabel('I'); ylabel('J'); zlabel('K')


%% From IMAGE:

% dim1 = I
% dim2 = L
% dim3 = A

% ijk2ras = [-3 -1 2]

[-1 -3 2]

%% From HEADER (transform matrix)

mni.vox2ras1
%     -2     0     0    76
%      0     0     2  -111
%      0    -2     0    87
%      0     0     0     1

%%
T = mni.vox2ras1;
V = permute(mni.vol,[2 1 3 4]);
ras = [nan nan 10]; % position x in WORLD-coordinate


R = inv(T);
[rasdim2ijkdim,~] = find(R(:,1:3));
[ijkdim2rasdim,~] = find(T(:,1:3)); % ijk2ras

% V = mni.vol;
rasd = ras;
rasd(isnan(ras)) = 0;
ijk = R*[rasd ones(size(ras,1),1)]'; % ras2ijk
rasdim = find(~isnan(ras));
ijkdim = rasdim2ijkdim(rasdim);

% World coordinates of IJK AXES
Xi = []; axisnames = 'ras';
for i = 1:3
  G = ones(size(V,i),4);
  G(:,i) = 1:size(V,i);
  U = T*G';
  Xi(i).axis = U(ijkdim2rasdim(i),:);
  Xi(i).axisname = axisnames(ijkdim2rasdim(i));
end

switch (ijkdim)
  case 1
    % Image is YX... not XY
    Vi = squeeze(interpn(V, ijk(idim), 1:size(V,2), 1:size(V,3)))';
    u = Xi(2);
    w = Xi(3);
  case 2
    Vi = squeeze(interpn(V, 1:size(V,1), ijk(idim), 1:size(V,3)))';
    u = Xi(1);
    w = Xi(3);
  case 3
    Vi = squeeze(interpn(V, 1:size(V,1), 1:size(V,2), ijk(idim)))';
    u = Xi(1);
    w = Xi(2);
end

clf
imagesc(u.axis, w.axis, Vi)
axis xy
xlabel(u.axisname); ylabel(w.axisname)
%%

i = [10.5 nan nan]
Vi = interp3(V, i(1), 1:size(V,2), 1:size(V,3));
clf
imagesc(x(:,1,1),squeeze(y(1,:,1)),Vi)


%%
clear
mri = load_nifti('~/MNI152_T1_1mm.nii')
[Vi,u,w] = vol2slice(mri.vol, mri.vox2ras, [nan 0.5 nan],'linaer');
imagesc(u.axis,w.axis,Vi); axis xy;
% xlabel(u.axisname); ylabel(w.axisname);
% axis([-20 20 20 60])
xlabel('MNI-Y = 0.5 mm'); grid off;

%%
mri = load_nifti('~/MNI152_T1_1mm.nii')

thal = load_nifti('/usr/local/fsl/data/atlases/Thalamus/Thalamus-prob-1mm.nii.gz')'

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
  
  axespos(ax,i)
  [Qi,u,w] = vol2slice(thal.vol(:,:,:,1), thal.vox2ras, [nan y nan],'linear');
  Qi(~Qi) = nan;
  h = pcolor(u.axis, w.axis, Qi);
  h.LineStyle = 'none';
  colormap(gca,turbo)
  axis off image
  
end

%%
mri = load_nifti('~/MNI152_T1_1mm.nii')

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

mri = load_nifti('~/MNI152_T1_1mm.nii');
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
mri = load_nifti('~/MNI152_T1_1mm.nii');
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

%% "Glass brain"
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
