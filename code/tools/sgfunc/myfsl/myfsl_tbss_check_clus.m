function JOB = myfsl_tbss_check_clus(JOB)

load /scr/vatikan3/APConn/mat/info17_sorted.mat
subjID = fsss_subjID(subjID17);
JOB.dir_tbss = '/scr/vatikan3/APConn/tbss_n17';
%if ~isfield(JOB,'dir_figure'), JOB.dir_figure=pwd; else [~,~]=mkdir(JOB.dir_figure); end
JOB.dir_figure='/scr/vatikan3/APConn/tbss_n17/fig_clus';
JOB.dir_dwi = '/scr/vatikan3/APConn/dwi1p7/';

%% first on the skeleton
model_desc = 'FA_1+AP+age+sex+asia_AP';
name_clus = 'AP_tfce_corrp_tstat1_clusterid_thrs_0.95.nii';
name_clusname = 'AP_tfce_corrp_tstat1_cluster_thrs_0.95.mat';
load ([JOB.dir_tbss,'/GLM/',model_desc,'/',name_clusname]);
clusname= clus;

if ~JOB.skip1
name_pval = 'AP_tfce_corrp_tstat1.nii.gz';
nii = load_uns_nii([JOB.dir_tbss,'/GLM/',model_desc,'/',name_clus]);
clus = nii.img;
nii = load_uns_nii([JOB.dir_tbss,'/GLM/',model_desc,'/',name_pval]);
pval = nii.img;
gunzip([JOB.dir_tbss,'/GLM/',model_desc,'/all_FA_skeletonised.nii.gz']);
numclus=max(clus(:));
for ci=1:numclus
hf = figure('position',[1923           2        1000        1150]);
c=clusname(ci).id;
ind = clus == c;
IJK = find3(ind);
ijk = mode(IJK);
ax = axeslayout1(20,[4 5],[.035 .07]);
meanFA=zeros(1,17);
FAatPeak=zeros(1,17);

i=1;
axespos(ax,i);
fa = load_uns_nii([JOB.dir_tbss,'/GLM/',...
model_desc,'/all_FA_skeletonised.nii'],i,[],[],[],[],ijk(3));
fa = fa.img;
favals =fa(ind(:,:,ijk(3)));
clim = [min(favals(:)), max(favals(:))];
imagesc(fa'); set(gca,'xdir','rev','ydir','nor'); 
colormap hot; axis image; axis off;
hold on;
contour(ind(:,:,ijk(3))',1,'color',[0 1 0],'linewidth',2);
title(['mnicoord=[',num2str(ijk2xyz(ijk,nii)),']mm']);
ax1=axes('position', [0.006 0.78 0.188 0.03]);
hb=colorbar('peer',ax1,'location','north');  axis(ax1,'off')
set(hb,'xcolor','w')
caxis(clim);

for i=1:17
axespos(ax,i+1);
fa = load_uns_nii([JOB.dir_tbss,'/GLM/',...
model_desc,'/all_FA_skeletonised.nii'],i,[],[],[],[],ijk(3));
fa = fa.img;
imagesc(fa'); set(gca,'xdir','rev','ydir','nor'); 
colormap hot; axis image; axis off;
hold on;
contour(ind(:,:,ijk(3))',1,'color',[0 1 0],'linewidth',1);
xmargin=10;
ymargin=20;
axis([min(IJK(:,1))-xmargin,max(IJK(:,1))+xmargin, ...
min(IJK(:,2))-ymargin, max(IJK(:,2))+ymargin]);
caxis(clim);
meanFA(i) = mean(fa(ind(:,:,ijk(3))));
maxfa= fa( pval(:,:,ijk(3))==max(max(pval(:,:,ijk(3)))));
FAatPeak(i) = max(maxfa);
title([subjID{i},':meanFA=',num2str(meanFA(i),2)]);
end

M = 1 + term(age17') + term(isFemale17') + term(isAsian17');
subplot(7,5,34)
SurfStatPlot(isAP17', meanFA'); xlim([-.5 1.5]); 
set(gca,'xtick',[0 1],'xticklabel',{'NAP','AP'}); ylabel('unadj. meanFA');

subplot(7,5,35)
SurfStatPlot(isAP17', FAatPeak', M); xlim([-.5 1.5]); 
set(gca,'xtick',[0 1],'xticklabel',{'NAP','AP'}); ylabel('adj. FA at Peak');

i=19;
axespos(ax,i);
caxis(clim);
text(0,0.1,{[clusname(ci).cog_name],...
['#vox=',num2str(sum(~~ind(:)))],...
['minP=',num2str(1-max(pval(ind)))]},'fontsize',15);
axis off

screen2png([JOB.dir_figure,'/tbss_cluster_check_skeleton_',...
string2filename(clusname(ci).cog_name),'.png']);
close(hf)
end
end

%% deprojected

if ~JOB.skip2
numclus=numel(clusname);
for ci=1:numclus
hf = figure('position',[1923           2        1000        1150]);
ax = axeslayout1(20,[4 5],[.035 .07]);
meanFA=zeros(1,17);
FAatPeak=zeros(1,17);

i=1;
fname_clus=[JOB.dir_tbss,'/FA/',subjID{i},'_FA_clusidmap.nii.gz'];
fname_fa  =[JOB.dir_tbss,'/FA/',subjID{i},'_FA.nii.gz'];
axespos(ax,i);
c=clusname(ci).id;
[clim,~,ijk] = showslice(fname_clus, fname_fa, c, 0);
%title(['mnicoord=[',num2str(ijk2xyz(ijk,nii)),']mm']);
ax1=axes('position', [0.006 0.78 0.188 0.03]);
hb=colorbar('peer',ax1,'location','north');  axis(ax1,'off')
set(hb,'xcolor','w')
clim=[clim(1)*0.9 clim(2)*1.1];
caxis(clim);

for i=1:17
axespos(ax,i+1);
fname_clus=[JOB.dir_tbss,'/FA/',subjID{i},'_FA_clusidmap.nii.gz'];
fname_fa  =[JOB.dir_tbss,'/FA/',subjID{i},'_FA.nii.gz'];
c=clusname(ci).id;
[~,meanFA(i)]=showslice(fname_clus, fname_fa, c, 1);
caxis(clim);
title([subjID{i},':meanFA=',num2str(meanFA(i),2)]);
end

M = 1 + term(age17') + term(isFemale17') + term(isAsian17');
subplot(7,5,34)
SurfStatPlot(isAP17', meanFA'); xlim([-.5 1.5]); 
set(gca,'xtick',[0 1],'xticklabel',{'NAP','AP'}); ylabel('unadj. meanFA');

subplot(7,5,35)
SurfStatPlot(isAP17', meanFA', M); xlim([-.5 1.5]); 
set(gca,'xtick',[0 1],'xticklabel',{'NAP','AP'}); ylabel('adj. meanFA');

i=19;
axespos(ax,i);
caxis(clim);
text(0,0.1,{[clusname(ci).cog_name]},'fontsize',15);
axis off

screen2png([JOB.dir_figure,...
'/tbss_cluster_check_deproj_',string2filename(clusname(ci).cog_name),'_fa.png']);
close(hf)
end
end


%% now vector map!
numclus=numel(clusname);
for ci=1:numclus
hf = figure('position',[1923           2        1000        1150]);
ax = axeslayout1(20,[4 5],[.035 .07]);
meanvec = zeros(17,3);

i=1;
fname_clus=[JOB.dir_dwi,'/',subjID{i},'.probtrkX/tbssAP/label_nn.nii.gz'];
fname_fa=[JOB.dir_dwi,'/',subjID{i},'/spm12_1p7_Tensor_FA.nii.gz'];
fname_vec =[JOB.dir_dwi,'/',subjID{i},'/spm12_1p7_Tensor_V1.nii.gz'];

axespos(ax,i);
c=clusname(ci).id;
[clim,~,ijk] = showslicevec1(fname_clus, fname_fa, fname_vec, c, 0);
title(['mnicoord=[',num2str(ijk2xyz(ijk,fname_vec)),']mm']);

for i=1:17
axespos(ax,i+1);
fname_clus=[JOB.dir_dwi,'/',subjID{i},'.probtrkX/tbssAP/label_nn.nii.gz'];
fname_fa=[JOB.dir_dwi,'/',subjID{i},'/spm12_1p7_Tensor_FA.nii.gz'];
fname_vec =[JOB.dir_dwi,'/',subjID{i},'/spm12_1p7_Tensor_V1.nii.gz'];
c=clusname(ci).id;
[~,meanvec(i,:)]=showslicevec1(fname_clus, fname_fa, fname_vec, c, 1);
title(sprintf('%s:v1=(%0.1f,%0.1f,%0.1f)', subjID{i}, ...
meanvec(i,1), meanvec(i,2), meanvec(i,3)));
end

for j=1:3
subplot(7,8,56-3+j)
[~,~,pval]=SurfStatPlot(isAP17', meanvec(:,j)); 
title(['p=',num2str(pval,1)])
xlim([-.5 1.5]); 
set(gca,'xtick',[0 1],'xticklabel',{'NAP','AP'}); ylabel(['unadj. V1',num2str(j)]);
end

i=19;
axespos(ax,i);
text(0,0.1,{[clusname(ci).cog_name]},'fontsize',15);
axis off

screen2png([JOB.dir_figure,'/tbss_cluster_check_deproj_',string2filename(clusname(ci).cog_name),'_v1.png']);
close(hf)
end



end


function [clim ,meanfa, ijk] = showslice(fname_clus, fname_fa, c, zoom, pval)


nii = load_uns_nii(fname_clus);
clus = nii.img;

ind = clus == c;
IJK = find3(ind);
ijk = mode(IJK);

fa = load_uns_nii(fname_fa,[],[],[],[],[],ijk(3));
fa = fa.img;
favals =fa(ind(:,:,ijk(3)));
clim = [min(favals(:)), max(favals(:))];
imagesc(fa'); set(gca,'xdir','rev','ydir','nor'); colormap hot; axis image; axis off;
hold on;
contour(ind(:,:,ijk(3))',1,'color',[0 1 0],'linewidth',2);
meanfa = max(favals);
% maxfa= fa( pval(:,:,ijk(3))==max(max(pval(:,:,ijk(3)))));
% maxfa = max(maxfa);

if zoom
xmargin=10;
ymargin=20;
axis([min(IJK(:,1))-xmargin,max(IJK(:,1))+xmargin, min(IJK(:,2))-ymargin, max(IJK(:,2))+ymargin]);
end

end

%%
function [clim ,meanvec, ijk] = showslicevec1(fname_clus, fname_fa, fname_v1, c, zoom)
clim=[];

nii = load_uns_nii(fname_clus);
clus = nii.img;

ind = clus == c;
IJK = find3(ind);
ijk = mode(IJK,1);

% DIFFERENT ORIENTATION!!!!!!!!!!!!
fa = load_uns_nii(fname_fa);%,[],[],[],[],[],ijk(3));
fa = squeeze(fa.img(:,ijk(2),:));
fa = zeroone(fa);

v1 = load_uns_nii(fname_v1);%,[],[],[],[],[],ijk(3));
v1 = squeeze(abs(squeeze(v1.img(:,ijk(2),:,:))));

favals =fa(ind(:,ijk(2),:));
%clim = [min(favals(:)), max(favals(:))];
slice=v1;
J=[3 1 2];
for j=1:3
slice(:,:,j) = v1(:,:,J(j)) .* fa;
end
imagesc(slice); axis image

%set(gca,'xdir','rev','ydir','nor'); colormap hot; axis image; axis off;
hold on;
contour(squeeze(ind(:,ijk(2),:)),1,'color',[1 1 1],'linewidth',2);
meanvec=[0 0 0];
for j=1:3
v1j = squeeze(v1(:,:,J(j)));
meanvec(j)=mean(v1j(ind(:,ijk(2),:)));
end
meanvec = meanvec / l2norm(meanvec);
% maxfa= fa( pval(:,:,ijk(3))==max(max(pval(:,:,ijk(3)))));
% maxfa = max(maxfa);

if zoom
xmargin=10;
ymargin=20;
axis([min(IJK(:,3))-xmargin,max(IJK(:,3))+xmargin, min(IJK(:,1))-ymargin, max(IJK(:,1))+ymargin]);
end

end


