function JOB = myspm_results_montage(JOB)
% JOB = myspm_results_montage(JOB)
%
% It creates summary montage for a given GLM directory
%
% JOB requires:
% (.dir_glm)    'Nx1' directory to save SPM results (default=pwd)
% (.slicedim)   [1xN] slice dimension (default=3, axial in RAS)
% (.fname_t1w)
% (.pitch)
% (.tight)      [1x1] crop images tightly to overlay (default=1)
%
% (cc) 2018, sgKIM. solleo@gmail.com, https://ggooo.wordpress.com/

if nargin < 1,  JOB=[];  end
if ~isfield(JOB,'dir_glm'), JOB.dir_glm=pwd; end

load([JOB.dir_glm,'/SPM.mat']);
% See if the even contrasts are flipped ones:
if sum(SPM.xCon(1).c + SPM.xCon(2).c) == 0
  bothsigns=1;
  K=1:2:numel(SPM.xCon);
else
  bothsigns=0;
  K=1:numel(SPM.xCon);
end
[~,modelname]=myfileparts(JOB.dir_glm);
JOB.fnames_png={};
untilK=numel(K);
if isfield(JOB,'untilK'), untilK=JOB.untilK; end
for i=1:untilK
  k=K(i);
  cfg=JOB;
  cfg.SPM=[SPM.xCon(k).STAT,'_',zeropad(k,4)];
  if bothsigns
    cfg.cntrst={SPM.xCon(k).name,SPM.xCon(k+1).name};
  else
    cfg.cntrst={SPM.xCon(k).name,''};
  end
  cfg.title=[modelname,':',SPM.xCon(k).name];
  if isfield(JOB,'title_suffix')
    cfg.title=[cfg.title,' ',JOB.title_suffix];
  end
  if isfield(cfg,'fname_struct') && ~isfield(cfg,'fname_t1w')
    cfg.fname_t1w=cfg.fname_struct;
  end
  if ~isfield(cfg,'tight'), cfg.tight=1; end
  JOB.fnames_png{i}=imagemontage_20181025(cfg);
end
end

function [fname,cfg]=imagemontage_20181025(cfg)
% cfg requires:
%  .dir_glm
% (.cntrst)  {1x2} dir_glm/sigclus_"+var".nii
%  .SPM       dir_glm/spm"T_0001".nii
%  .title
% (.mni_xyz)
%  .slicedim
% (.radius)
% (.seedcoord)
% (.fname_seed)
% (.fname_t1w)
% (.baseimagerange)
%
% (cc) 2016, sgKIM

fname='';
if ~nargin, help(mfilename),return;end
[~,pathname]=myfileparts(cfg.dir_glm);
if ~isfield(cfg,'radius'), radius=10; else radius=cfg.radius; end
pwd0=pwd;
cd(cfg.dir_glm);
if ~isfield(cfg,'dir_fig')
  dir_fig=pwd;
else
  dir_fig=cfg.dir_fig;
end
if ~isfield(cfg,'fname_t1w')
  fname_t1w='~/MNI152_T1_1mm.nii';
else
  fname_t1w=cfg.fname_t1w;
end
cmap=sgcolormap('GRAY-BLUE-RED-BRIGHT');
if isfield(cfg,'colormap')
  cmap=colormap;
end
t1w=load_uns_nii(fname_t1w);
if strcmp(fname_t1w,'~/MNI152_T1_1mm.nii')
  baseimagerange=[2500 11000];
else
  imrange=[mode(t1w.img(:)) max(t1w.img(:))];
  baseimagerange=[imrange(1)+diff(imrange)*0.1 imrange(2)];
end
if isfield(cfg,'baseimagerange')
  baseimagerange=cfg.baseimagerange;
end

%% Read T(F)-map and sigclus-map(s)
if isfield(cfg,'custommap')
  fname_nii=cfg.custommap;
else
  fname_nii=['spm',cfg.SPM,'.nii'];
end
if ~isfield(cfg,'cntrst')
  cfg.cntrst={'',''};
end
over=load_untouch_nii_like(fname_nii,fname_t1w);
fname1=['sigclus_',cfg.cntrst{1},'.nii'];
fname2=['sigclus_',cfg.cntrst{2},'.nii'];
if exist(fname1,'file')
  hdr=load_untouch_header_only(fname1);
  clus1=load_untouch_nii_like(fname1, fname_t1w); % interpolation for visualization
else
  clus1=over;
  clus1.img=clus1.img*0;
end
if exist(fname2,'file')
  hdr=load_untouch_header_only(fname2);
  clus2=load_untouch_nii_like(fname2, fname_t1w); % interpolation for visualization
else
  clus2=over;
  clus2.img=clus2.img*0;
end
if exist(fname1,'file') || exist(fname2,'file')
  thres_u = str2double(hdr.hist.descrip(22:26));
  thres_k = str2double(hdr.hist.descrip(33:end));
  if isfield(cfg,'caxis')
    caxis1=cfg.caxis;
  else
    caxis1=[thres_u max(abs(over.img(:)))];
  end
  clus=clus1.img>0.5;
  clus(clus2.img>0.5) = 1;
  clus=double(clus);
  T_thres = over.img .* ~~clus;
else
  warning('No significant clusters/voxels found. Nothing to display.')
  return
end
%% Read a table
tab = myspm_readtable;
if ~isempty(tab)
  idx=ismember(tab.contrast_name, cfg.cntrst{1});
  mni_xyz=[tab.x_mm(idx) tab.y_mm(idx) tab.z_mm(idx)];
  idx=ismember(tab.contrast_name, cfg.cntrst{2});
  mni_xyz=[mni_xyz; tab.x_mm(idx) tab.y_mm(idx) tab.z_mm(idx)];
end
if isfield(cfg,'mni_xyz'), mni_xyz=cfg.mni_xyz; end
mni_ijk=round(xyz2ijk(mni_xyz,over));
if ~isfield(cfg,'slicedim')
  slicedim=ones(1,size(mni_ijk,1))*3;
else
  slicedim=cfg.slicedim;
end
tab_mni_ijk=mni_ijk(1,:);
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
if numel(slicedim) == 1
  slicedim=repmat(slicedim,[1 num_clus]);
end
base=double(t1w.img);
overlay=double(over.img);
%% rotation images?
if isfield(cfg,'pitch')&&cfg.pitch
  base = imrotate3(base,cfg.pitch,[0 1 0],'cubic','crop','FillValues',0);
  T_thres = imrotate3(T_thres,cfg.pitch,[0 1 0],'cubic','crop','FillValues',0);
  overlay = imrotate3(overlay,cfg.pitch,[0 1 0],'cubic','crop','FillValues',0);
end
%% Bounding box
IJK=find3(~~overlay);
bb=[min(IJK)-1; max(IJK)+1];

%% main slices
if isfield(cfg,'tight')
  FOV=abs(diff(bb,1));
else
  FOV=size(overlay);
end

ax=axeslayout_(num_clus+1,[1 slicedim],FOV);
figpos=[-1721 127 sum(ax.w0)/FOV(1)*120 147/0.85];
figure('color',[0 0 0], 'position',figpos);
for i=1:num_clus
  switch slicedim(i)
    case 3
      cfg1 = struct('ijk',mni_ijk(i,:));
      label1=['z = ',num2str(mni_xyz(i,3)),' mm'];
    case 2
      cfg1 = struct('ijk',mni_ijk(i,:));
      label1=['y = ',num2str(mni_xyz(i,2)),' mm'];
    case 1
      cfg1 = struct('ijk',mni_ijk(i,:));
      label1=['x = ',num2str(mni_xyz(i,1)),' mm'];
  end
  j=i+1;
  axespos(ax,j);
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
  if isfield(cfg,'subcort')
    cfg1.subcort=cfg.subcort(i);
  end
  if isfield(cfg,'subcortbox')
    cfg1.subcortbox=cfg.subcortbox(i);
  end
  if isfield(cfg,'subcortzoom')
    cfg1.subcortzoom=cfg.subcortzoom(i);
  end
  if isfield(cfg,'subcortannot')
    cfg1.subcortannot = cfg.subcortannot(i);
  end
  if isfield(cfg,'roicontourcolor')
    cfg1.roicontourcolor = cfg.roicontourcolor;
  end
  if isfield(cfg,'tight')
    cfg1.bb = bb;
  end
  cfg1.colormap=cmap;
  
  % execute imageover1
  imageover1(base, T_thres, cfg1);
  
  if slicedim(i)<3
    axis0=axis;
    axis0(3)=axis0(3)-(FOV(2)-FOV(3));
    axis(axis0)
  end
  
  % coordinate label
  hal=axespos(ax,j);
  if ismac
  coordinate_label_fontsize=10;
  elseif isunix
      coordinate_label_fontsize=7;
  end
  text(hal,0.5,0.05,label1,'fontsize',coordinate_label_fontsize, ...
    'backgroundcolor','none', 'color','w', ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle');
  axis off
end

%% Saggittal navigator
axespos(ax,1)
if ~isfield(cfg,'sagnav_xyz')
  cfg1 = struct('ijk',tab_mni_ijk, 'slicedim',1, 'caxis',caxis1);
else
  cfg1 = struct('ijk',round(xyz2ijk(cfg.sagnav_xyz,t1w)), ...
    'slicedim',1, 'caxis',caxis1);
end
N=size(mni_ijk,1);
if isfield(cfg,'lines_y')
  cfg1.lines_y = cfg.lines_y;
else
  cfg1.lines_y = slicedim==2;
end
for j=1:N
  if slicedim(j)==2 && cfg1.lines_y(j)
    cfg1.slices(j) = mni_ijk(j,2);
  elseif slicedim(j)==3
    cfg1.slices(j) = mni_ijk(j,3);
  else
    cfg1.slices(j) = nan;
  end
end
if isfield(cfg,'fname_seed') && ~isfield(cfg,'roi')
  cfg1.roi=roi.img; % seed mask contour
elseif isfield(cfg,'roi')
  cfg1.roi=double(roi);
end
cfg1.baseimagerange = baseimagerange;
if isfield(cfg,'baseedge')
  cfg1.baseedge = cfg.baseedge;
end
if isfield(cfg,'nolines'), cfg1.nolines=cfg.nolines; end
cfg1.colormap=cmap;
imageover1(base, T_thres, cfg1);
axis0=axis;
axis0(3)=axis0(3)-(FOV(2)-FOV(3));
axis(axis0)

% % coordinate label
% hal=axespos(ax,j);
% xyz_sag=round(ijk2xyz(cfg1.ijk, t1w));
% label1=['x = ',num2str(xyz_sag(1,1)),' mm'];
% coordinate_label_fontsize=8;
% text(hal,0.5,0.05,label1,'fontsize',coordinate_label_fontsize, ...
%   'backgroundcolor','none', 'color','w', ...
%   'HorizontalAlignment','center', 'VerticalAlignment','middle');
% axis off

if ~isfield(cfg,'notitle')
  if ~isfield(cfg,'title')
    cfg.title=pathname;
  end
  Title=cfg.title;
  if ~isfield(cfg,'titlefontsize')
    if ismac
    cfg.titlefontsize=19; 
    elseif isunix
      cfg.titlefontsize=14;
    end
  end
  % TITLE label
  hal=axes('position',[0 0 1 1]); axis off;
  text(hal,0.003,0.93, Title, 'fontsize',cfg.titlefontsize, ...
    'color','w', 'interp','none', 'BackgroundColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment','Middle');
end

%% colorbar
ax2=ax;
ax2.h(1)=(FOV(2)-FOV(3))/FOV(2)*ax2.h(1);
ax2.w(1)=ax2.w(1)/2;
ax2.x(1)=ax2.w(1);
H=axespos(ax2,1);
set(H,'CLim',[caxis1(1) caxis1(2)],'color','k');
xtick=[caxis1(1)+diff(caxis1)*0.1 caxis1(2)-diff(caxis1)*0.1];
xticklabel={sprintf('%.1f',caxis1(1)), sprintf('%.1f',caxis1(2))};
cb1=colorbar('peer',H,'north','color','w','xtick',xtick,'xticklabel',xticklabel);
colormap(H,cmap(129:end,:)) % red-yellow
pos=get(cb1,'position');
set(cb1,'position',[pos(1) pos(2)+pos(4)*1.5 pos(3) pos(4)])

ax2=ax;
ax2.h(1)=(FOV(2)-FOV(3))/FOV(2)*ax2.h(1);
ax2.w(1)=ax2.w(1)/2;
H=axespos(ax2,1);
caxis2=[-caxis1(2) -caxis1(1)];
set(H,'CLim',[caxis2(1) caxis2(2)],'color','k');
xtick=[caxis2(1)+diff(caxis2)*0.15 caxis2(2)-diff(caxis2)*0.15];
xticklabel={sprintf('%.1f',caxis2(2)), sprintf('%.1f',caxis2(1))};
cb2=colorbar('peer',H,'north','color','w','xdir','rev',...
  'xtick',xtick,'xticklabel',xticklabel);
colormap(H,cmap(65:128,:)) % blue-cyan
pos=get(cb2,'position');
set(cb2,'position',[pos(1) pos(2)+pos(4)*1.5 pos(3) pos(4)])

%% save
if isfield(cfg,'suffix')
  suffix=cfg.suffix;
else
  suffix='';
end
if ~isfield(cfg,'fname_png')
  fname=[dir_fig,'/fig_',pathname,'_',cfg.cntrst{1},suffix,'.png'];
else
  fname=cfg.fname_png;
end
[d1,~,~]=fileparts(fname);
[~,~]=mkdir(d1);
if ~isfield(cfg,'dpi')
  dpi=150;
else
  dpi=cfg.dpi;
end
screen2png(fname,dpi,1);

cd(pwd0);
end

function axes_pos = axeslayout_(N,slicedim,ijkdim)
h=0.85*ones(1,N);
y=zeros(1,N);
w_i=slicedim*0;
w_i(slicedim==1)=2;
w_i(slicedim==2)=1;
w_i(slicedim==3)=1;
w0=ijkdim(w_i);
w=w0./sum(w0);
x(1)=0;
for i=2:N
  x(i)=x(i-1)+w(i-1);
end
axes_pos = struct('x',x,'y',y,'w',w,'h',h,'w0',w0);
end
