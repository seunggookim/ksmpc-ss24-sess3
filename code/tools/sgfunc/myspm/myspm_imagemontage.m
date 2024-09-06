function fname=myspm_imagemontage(cfg)
% cfg requires:
% .pathname,
% .cnt
% .mni_xyz
% .slicedim
% .radius
% .seedcoord
%
% (cc) 2016, sgKIM

if ~nargin, help(mfilename); return; end

dir_glm = cfg.dir_glm;
if isfield(cfg,'mni_xyz'), mni_xyz = cfg.mni_xyz; end

% optional
if ~isfield(cfg,'radius'), radius=10; else radius=cfg.radius; end
if ~isfield(cfg,'slicedim'), slicedim=3; else slicedim=cfg.slicedim; end
if ~isfield(cfg,'seedcoord'), nocrosshair=1; seedcoord=0; else seedcoord=1; end

caxis1=[-7 7];
caxis2=[-4 4];
if isfield(cfg,'threshold');
caxis1=[-cfg.threshold(2) cfg.threshold(2)];
caxis2=[-cfg.threshold(1) cfg.threshold(1)];
end

pwd0=pwd;
load /scr/vatikan4/conmus3/mat/info23.mat
if ~isfield(cfg,'dir_fig')
dir_fig=dir_glm;
else
dir_fig=cfg.dir_fig;
end
baseimagerange=[2500 10500];
if ~isfield(cfg,'fname_t1w')
fname_t1w='/usr/share/fsl/5.0/data/standard/MNI152lin_T1_1mm_brain.nii.gz';
else
fname_t1w=cfg.fname_t1w;
end
% 70% intensity for z=97:end
t1w=load_uns_nii(fname_t1w);
if ~isfield(cfg,'SPM'), cfg.SPM='T'; end
switch cfg.SPM
case 'T'
[~,fnameT]=mydir([dir_glm,'/spmT_0001.*']);
over =load_untouch_nii_like(fnameT,fname_t1w);
case 'con'
[~,fnameT]=mydir([dir_glm,'/con_0001.*']);
over =load_untouch_nii_like(fnameT,fname_t1w);
end
if ~isfield(cfg,'var')
fname1=[dir_glm,'/sigclus_+1.nii'];
else
fname1=[dir_glm,'/sigclus_+',cfg.var,'.nii'];
end
if exist(fname1,'file')
clus1=load_untouch_nii_like(fname1, fname_t1w);
% interpolation for visualization! (nice)
else
clus1=over; clus1.img=clus1.img*0;
end
if ~isfield(cfg,'var')
fname2=[dir_glm,'/sigclus_-1.nii'];
else
fname2=[dir_glm,'/sigclus_-',cfg.var,'.nii'];
end
if exist(fname2,'file')
clus2=load_untouch_nii_like(fname2, fname_t1w);
else
clus2=over; clus2.img=clus2.img*0;
end
if ~exist(fname1,'file') && ~exist(fname2,'file')
warning('No significant cluster found')
fname=[];
return
end
clus=clus1.img;
clus(~~clus2.img) = clus2.img(~~clus2.img) + max(clus1.img(:));
clus=round(clus);

T_thres = over.img .* ~~clus;
cd(dir_glm);
tab = myspm_readtable;
if ~exist('mni_xyz','var')||isempty(mni_xyz)
mni_xyz=tab.mni_xyz;
if slicedim == 3
z=unique(mni_xyz(:,3));
mni_xyz=[zeros(numel(z),2) z];
elseif slicedim == 2
y=unique(mni_xyz(:,2));
mni_xyz=[zeros(numel(y),1) y zeros(numel(y),1)];
end
end
tab_mni_ijk = xyz2ijk(tab.mni_xyz(1,:),over);
mni_ijk=xyz2ijk(mni_xyz,over);
if numel(mni_ijk) == 0; return; end
if ~isfield(cfg,'seedcoord')
seedcoord=mni_xyz(1,:);
else
seedcoord=cfg.seedcoord;
end
seedcoord_ijk = xyz2ijk(seedcoord, over);
num_clus=size(mni_xyz,1);

%% seed ROI
if isfield(cfg,'fname_seed') && ~isfield(cfg,'roi')
roi = load_untouch_nii_like(cfg.fname_seed, fname_t1w);
roi.img = round(roi.img);
elseif isfield(cfg,'roi')
roi = cfg.roi;
end

%%
if slicedim(1) <= 2
figpos=[54   318   1020   144]; % for six slices...
elseif slicedim(1) == 3
figpos=[54   318   968   201];
end
if isfield(cfg,'width_b')
figpos(3)=cfg.width_b;
end

figure('color','k', 'position',figpos);

%% Saggittal navigator
ax=axeslayout1(num_clus+1,[1,num_clus+1]);
axespos(ax,1)
cfg1 = struct('ijk',tab_mni_ijk, 'slicedim',1, 'caxis',caxis1); % sagittal section
% 2016-12-25: coronal section markers
N=size(mni_ijk,1);
if ~isfield(cfg,'lines_y')
cfg1.lines_y=zeros(1,N);
elseif numel(cfg.lines_y) == 1
cfg1.lines_y=ones(1,N) * cfg.lines_y;
else
cfg1.lines_y=cfg.lines_y;
end
for j=1:N
if cfg1.lines_y(j)
cfg1.slices(j) = mni_ijk(j,2);
else
cfg1.slices(j) = mni_ijk(j,3);
end
end
cfg1.slices
cfg1.baseimagerange=baseimagerange;
if isfield(cfg,'baseedge')
cfg1.baseedge = cfg.baseedge;
end
imageover1(double(t1w.img), T_thres, cfg1);
axis1=axis;
if ~isfield(cfg,'notitle')
cnt=cfg.cnt;
if ~exist('cnt','var')
cnt=dir_glm;
end
if ~iscell(cnt) && strcmp(cnt(3),'4')
cnt=[cnt(1:2),'-',cnt(5:6)];
end
text(mean(axis1(1:2)), axis1(4)+60, cnt, 'fontsize',14, 'color','w', ...
'HorizontalAlignment', 'center', 'VerticalAlignment','top', ...
'interp','none', 'BackgroundColor', 'k');
if isfield(cfg,'panidx')
text(axis1(2)-5,axis1(4)+45, cfg.panidx,'color','w','fontsize',16);
end
end
if numel(slicedim) == 1
slicedim=repmat(slicedim,[1 num_clus]);
end

%% main slices
for i=1:num_clus
axespos(ax,i+1)
if slicedim(i) == 3
cfg1 = struct('ijk',mni_ijk(i,:), 'label',['z = ',num2str(mni_xyz(i,3)),' mm']);
elseif slicedim(i) == 2
cfg1 = struct('ijk',mni_ijk(i,:), 'label',['y = ',num2str(mni_xyz(i,2)),' mm']);
end
if isfield(cfg,'seedcoord')
cfg1.crosshair=seedcoord_ijk;
cfg1.sphere.center_vox = seedcoord_ijk;
cfg1.sphere.radius_vox = radius/mean(t1w.hdr.dime.pixdim(2:4));
end
cfg1.xdir='rev';
cfg1.baseimagerange=baseimagerange;
cfg1.caxis=caxis1;
cfg1.slicedim = slicedim(i);
if i == 1
cfg1.isR=1;
end
if isfield(cfg,'isR')
cfg1.isR=cfg.isR(i);
end
if isfield(cfg,'fname_seed') && ~isfield(cfg,'roi')
cfg1.roi=roi.img; % seed mask contour
elseif isfield(cfg,'roi')
cfg1.roi=double(roi);
end
if isfield(cfg,'baseedge')
cfg1.baseedge = cfg.baseedge;
end
imageover1(double(t1w.img), T_thres, cfg1);
end
if num_clus <=5
xshift=0.025; yshift=-0.06;
else
xshift=0.03; yshift=-0.06;
end

%% colorbar
if strcmp(version,'8.2.0.701 (R2013b)')
ax2=ax;
ax2.y(1)=0.8;
ax2.h=0.2;
pos1=[ax2.w*0.1    0.8557-0.13    ax2.w*0.8   0.0311];
pos2=[ax2.w*0.1    0.7300-0.13    ax2.w*0.8   0.0311];
H=axespos(ax2,1); hold on;
set(H,'color','k', 'visible','off')

cb1=colorbar('peer',H,'south','color','w');
set(cb1,'position',pos1, 'color','w')
cmap=sgcolormap('GRAY-BLUE-RED-LIGHT');
if isfield(cfg,'colormap')
cmap=colormap;
end
caxis([caxis2(2) caxis1(2)]);
xticks=get(cb1,'xtick');
if sum(T_thres(:)>0)
colormap(cmap(129:end,:)) % red-yellow
set(cb1, 'xticklabel', xticks);
elseif sum(-T_thres(:)>0)
colormap(cmap(65:128,:)) % blue-cyan
set(cb1, 'xticklabel', -xticks);
end
cbfreeze(cb1)

if sum(-T_thres(:)>0) && sum(T_thres(:)>0)
cb2=colorbar('peer',H,'south','color','w');
caxis([-caxis2(1) -caxis1(1)]);
colormap(cmap(65:128,:)) % blue-cyan
set(cb2,'position',pos2, 'color','w')
xticks=get(cb2,'xtick');
set(cb2, 'xticklabel', -xticks);
cbfreeze(cb2)
end

axes(H)
if isfield(cfg,'mapname')
text(0.5, 0.5, cfg.mapname, 'fontsize',15, 'color','w', ...
'HorizontalAlignment', 'center');
end
if isfield(cfg,'threstxt')
text(0.5, -3.15, cfg.threstxt, 'fontsize',10, 'color','w', ...
'HorizontalAlignment', 'center');
end
else
warning('matlab version 8.2.0.701 (R2013b) is required to create colorbars');
end

%% save
if isfield(cfg,'suffix'), suffix=cfg.suffix; else, suffix='';
end
if ~isfield(cfg,'fname_fig')
fname=[dir_fig,'/fig.png'];
else
fname=cfg.fname_fig;
end
[d1,~,~]=fileparts(fname);
[~,~]=mkdir(d1);
screen2png(fname,900);
cd(pwd0);
end
