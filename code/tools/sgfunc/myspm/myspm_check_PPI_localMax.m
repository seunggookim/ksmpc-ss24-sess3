function myspm_check_PPI_localMax(JOB)
% myspm_check_PPI_localMax(JOB)
%
% JOB
%  .dir_glm
%  .dir_ppi
% checks original coordinates and adjusted local maxima
%
% (cc) 2015, sgKIM.   solleo@gmail.com   https://ggooo.wordpress.com

dir0 = JOB.dir_base;
JOB.subjID = fsss_subjID(JOB.subjID);
NumSubj = numel(JOB.subjID);

figure('color','k', 'position',[1921, 1, 1024, 1176]);
[~,ax] = getLayout(NumSubj+1);
for i = 1 : NumSubj
subjid = JOB.subjID{i};  
fname_strc = fullfile(dir0,subjid,JOB.name_strc);
fname_func = fullfile(dir0,subjid,JOB.dir_glm,JOB.name_func);
hdr_strc = load_nii_hdr(fname_strc);
hdr_func = load_nii_hdr(fname_func);

if sum(hdr_strc.dime.dim(2:4) == hdr_func.dime.dim(2:4)) ~= 3 % when dimensions are not matched
[~,b,c] = fileparts(JOB.name_func);
fname_strc_in_func = [dir0,'/',subjid,'/',b,'_in_func',c];
if ~exist(fname_strc_in_func,'file')
unix(['mri_convert --like ',fname_func,' ',fname_strc,' ',fname_strc_in_func]);
end
fname_strc = fname_strc_in_func; % just resample functional images for visualization
end
base = load_nii(fname_strc);
func = load_nii(fname_func);
if ~isfield(JOB,'dir_ppi') % only for Physiological voi check
JOB.dir_ppi = JOB.dir_glm;
end
load(fullfile(dir0,subjid,JOB.dir_ppi,'PPI.mat'),'PPI');

d = size(base.img);
coord_xyz =  PPI.xY.xyz'; % read from PPI.mat
coord_ijk =  round(xyz2ijk(coord_xyz,fname_func));
coord_ijk0 = round(xyz2ijk(JOB.voi.coord,fname_func)); % what I entered

axespos(ax,i); hold on;
img1 = base.img(:,:,coord_ijk(3))';
img2 = func.img(:,:,coord_ijk(3))';
img2(isnan(img2))=0;
img2(img2<JOB.thres) = 0;
imagesc(zeroone(img1)+zeroone(img2));
cmap=sgcolormap('GRAY-RED');
if isfield(JOB,'isneg')&&JOB.isneg, cmap=sgcolormap('GRAY-BLUE'); end
colormap(cmap);
caxis([0 2])

% mark the original coordinate
line([d(1)-coord_ijk0(1);d(1)-coord_ijk0(1)],[0;d(2)+1],'color','w');
line([0;d(1)+1], [coord_ijk0(2);coord_ijk0(2)],'color','w');

% if different, mark the shifted coodinate
if sum(coord_ijk0==coord_ijk)<3, % if the coordinates are not identical
vx=d(1)-coord_ijk(1); vy=coord_ijk(2); vr=JOB.voi.radius;
rectangle('position',[vx-vr, vy-vr, vr*2, vr*2],'curvature',[1 1],'edgecolor','k');
end
axis([1 d(1) 1 d(2)]); axis off; axis image; set(gca,'ydir','nor');
if isfield(JOB,'masksubjid')&&JOB.masksubjid
subjname=['subj',pad(i,3)];
else
subjname=subjid;
end
text(23,10,[subjname,':T>',num2str(round(JOB.thres*10)/10)], 'backgroundColor','w')

drawnow;
end
axespos(ax,NumSubj+1);
set(gca,'color','k'); axis off;
text(0,0.3,JOB.voi.name,  'interp','none','color','w');
text(0,0.2,JOB.dir_glm,   'interp','none','color','w');
text(0,0.1,JOB.name_func, 'interp','none','color','w');
screen2png([JOB.dir_fig,'/',JOB.voi.name,'_',JOB.dir_glm,'_',JOB.name_func,'.png']);

close(gcf)
end
