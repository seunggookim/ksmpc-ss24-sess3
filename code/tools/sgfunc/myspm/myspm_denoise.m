function JOB = myspm_denoise(JOB)
% JOB = myspm_denoise(JOB)
%
% creates diagnostic figures from the results of myy_compcor.m
%
% JOB requires:
%  .dir_data
%  .name_epi;
%  .TR_sec
%  .t1w_suffix
% (.num_pcs)
% (.cov_idx)    1x1 for rp+cc+art+gs
% (.subjid) % for a figure...
% (.name_rp)
% (.nofigure)
%
% (cc) 2015, 2017, sgKIM.

if ~isfield(JOB,'fig_dpi'),fig_dpi=300;else,fig_dpi=JOB.fig_dpi;end
% [TODO] option for sagittal slice!
% [TDOO] topological measures? degree map? efficiency boxplot?

% 0. set parameters
if ~isfield(JOB,'sag'), JOB.sag=0; end
path0=pwd;
TR_sec = JOB.TR_sec;
sumstat=[];
if ~isfield(JOB,'num_pcs'), JOB.num_pcs=16; end
if ~isfield(JOB,'detrend'), JOB.detrend=1;  end
if ~isfield(JOB,'varnorm'), JOB.varnorm=1;  end
if ~isfield(JOB,'bpf'),    JOB.bpf=[0,Inf]; end
if isfield(JOB,'restbpf'), JOB.bpf=[0.009 0.08]; end
FilterBand = JOB.bpf;
if ~isfield(JOB,'param_cc')
 JOB.param_cc = sprintf('n%db%0.2f-%0.2f',JOB.num_pcs, JOB.bpf);
end

if ~isfield(JOB,'param_art')
 if isfield(JOB,'global_threshold') && isfield(JOB,'motion_threshold')
  JOB.param_art = sprintf('%0.1fstd_%0.1fmm', ...
   JOB.global_threshold, JOB.motion_threshold);
 else
  JOB.param_art='3.0std_0.5mm';
 end
end
[~, name1, ext1]=fileparts(JOB.name_epi);
if ~isfield(JOB,'name_rp')
 [~,name1,~]=fileparts(JOB.name_epi);
 JOB.name_rp=['rp_',name1,'.txt'];
end
if ~isfield(JOB,'name_cc')
 JOB.name_cc=[name1,'_',JOB.param_cc,'_eigenvec.txt'];
end
if ~isfield(JOB,'t1w_suffix'), JOB.t1w_suffix='Brain'; end
if ~isfield(JOB,'fname_gmmask')
 JOB.fname_gmmask=[name1,'_c1',JOB.t1w_suffix,'.nii'];
end

path1=[fullfile(JOB.dir_data),'/'];
cd(path1);
if ~isfield(JOB,'dir_func')
 %JOB.dir_func=fullfile(JOB.dir_data,'fig');
 JOB.dir_func=fullfile(JOB.dir_data);
end
[~,~]=mkdir(JOB.dir_func);

%% read data
JOB.fname_epi = fullfile(JOB.dir_data, JOB.name_epi);
nii = load_uns_nii(JOB.fname_epi);
d = nii.hdr.dime.dim(2:5);
y = double(reshape(nii.img,[],d(4))');

%% find compcor regressors
fname_txt=[path1,JOB.name_cc];
JOB.fname_out = fname_txt;
%if ~exist(fname_txt,'file')
myy_compcor(JOB);
%end
cc  = load(fname_txt);
%% find rp regressors
fname_rp=JOB.name_rp;
rp=load(fname_rp); % [mm rad];
rp=[rp, ...
 l2norm([0 0 0; diff(rp(:,1:3))]) , ...
 l2norm([0 0 0; diff(rp(:,4:6))]) ]; % 2019-06-11
[p6,f6,e6] = fileparts(fname_rp);
save([p6,f6,'8',e6],'rp','-ascii');

% find global signal
gs = load([path1,'/',name1,'_gm.txt']);
output_suffix=['_',JOB.param_cc];

% find the center of precuneous
if isfield(JOB,'dir_fs') && isfield(JOB,'subjid')
 subjid=JOB.subjid;
 fname_aparc=fullfile(JOB.dir_data,'oaparc.a2009s+aseg.nii');
 if ~exist(fname_aparc,'file')
  unix(['mri_convert --resample_type nearest ', ...
   fullfile(JOB.dir_fs,subjid,'mri','aparc.a2009s+aseg.mgz'), ...
   ' ', fname_aparc]);
 end
 aparc = load_uns_nii(fname_aparc);
 prec=[];
 aseglabels=[12130 11130]; % precuneous
 %aseglabels=[12109 12109]; % post cingulate gyrus
 prec.ijk = find3(ismember(aparc.img,aseglabels));
 prec.ijk_med = round(median(prec.ijk));
 
 if ~ismember(aparc.img(prec.ijk_med(1),prec.ijk_med(2),prec.ijk_med(3)),aseglabels)
  % if not, find nearest precuneous voxel
  prec.dist=l2norm(prec.ijk - repmat(prec.ijk_med,[size(prec.ijk,1),1]));
  [~,prec.idx] = min(prec.dist);
  prec.ijk_med = prec.ijk(prec.idx,:);
  % double check
  if ~ismember(aparc.img(prec.ijk_med(1),prec.ijk_med(2),prec.ijk_med(3)),aseglabels)
   error('??!!')
  end
 end
 
 prec.xyz_med = ijk2xyz(prec.ijk_med, aparc);
 prec.xyz = ijk2xyz(prec.ijk, aparc);
 clear aparc
 
 if JOB.sag
  % sagittal slice
  crrmap=[];
  % find 1-d index of precuneous mask
  crrmap.ijk = round(xyz2ijk(prec.xyz, nii)); % find voxel-coord of the precuneus
  crrmap.ind_prec  = sub2ind(d(1:3), crrmap.ijk(:,1), crrmap.ijk(:,2),crrmap.ijk(:,3)); % find 1-d index of precuneus
  crrmap.mask_prec = nii.img(:,:,:,1)*0; % 3-d volume
  crrmap.mask_prec(crrmap.ind_prec) = true;
  
  % find 1-d index of a slice at the median of precuneus
  crrmap.ijk_med = round(xyz2ijk(prec.xyz_med, nii)); % find voxel-coords of median of the precuneus
  crrmap.ind = sub2ind(d(2:3), crrmap.ijk_med(2), crrmap.ijk_med(3)); % find 1-d index of the median coord from the sagittal slice
  
  % find 1-d index of a slice
  crrmap.mask = nii.img(:,:,:,1)*0; % 3-d volume
  crrmap.mask(crrmap.ijk_med(1),:,:)=true; % find an axial slice at the median
  crrmap.idx = find(crrmap.mask(:)); % find 1-d index of the sagittal slice
  crrmap.y = y(:,crrmap.idx); % get GM values from voxels of the slice
  
  % find 1-d index of precuneus mask in the selected slice
  crrmap.mask_prec = crrmap.mask_prec & crrmap.mask; % precuneus AND an axial slice
  crrmap.idx_prec = find(crrmap.mask_prec(crrmap.ijk_med(1),:,:));
  
 else   % axial slice
  crrmap=[];
  % find 1-d index of precuneous mask
  crrmap.ijk = round(xyz2ijk(prec.xyz, nii)); % find voxel-coord of the precuneus
  crrmap.ind_prec  = sub2ind(d(1:3), crrmap.ijk(:,1), crrmap.ijk(:,2),crrmap.ijk(:,3)); % find 1-d index of precuneus
  crrmap.mask_prec = nii.img(:,:,:,1)*0; % 3-d volume
  crrmap.mask_prec(crrmap.ind_prec) = true;
  
  % find 1-d index of a slice at the median of precuneus
  crrmap.ijk_med = round(xyz2ijk(prec.xyz_med, nii)); % find voxel-coords of median of the precuneus
  crrmap.ind = sub2ind(d(1:2), crrmap.ijk_med(1), crrmap.ijk_med(2)); % find 1-d index of the median coord from the sagittal slice
  
  % find 1-d index of a slice
  crrmap.mask = nii.img(:,:,:,1)*0; % 3-d volume
  crrmap.mask(:,:,crrmap.ijk_med(3))=true; % find an axial slice at the median
  crrmap.idx = find(crrmap.mask(:)); % find 1-d index of the sagittal slice
  crrmap.y = y(:,crrmap.idx); % find GM values
  
  % find 1-d index of precuneus mask in the selected slice
  crrmap.mask_prec = crrmap.mask_prec & crrmap.mask; % precuneus AND an axial slice
  crrmap.idx_prec = find(crrmap.mask_prec(:,:,crrmap.ijk_med(3)));
 end
end

% get 256 GM voxels from the slice?
%gm = load_uns_nii(JOB.fname_gmmask);
gm.hdr  = load_untouch_header_only(JOB.fname_gmmask);
gm.img  = spm_read_vols(spm_vol_nifti(JOB.fname_gmmask));
gm.img(isnan(gm.img))=0;
idx = find(~~gm.img(:));
idx = idx(round(linspace(1,numel(idx),256)));
y=y(:,idx); % (regularly sampled) GM voxels
nii.img=[];


%% plot#1: timeseries
hf=figure('position',[1 48 958 1092],'visible','off');
set(gcf,'DefaultAxesColorOrder',get(0,'FactoryAxesColorOrder'))
ax0=subplot(6,7,[1:6]);
x1=[1:size(rp,1)];
plot(x1, [rp(:,1:3), rp(:,4:6)/pi*180])
set(gca,'xlim',[1,size(rp,1)]);
axespos = get(gca,'position');
legend({'x [mm]','y [mm]','z [mm]','rl [deg]','ph [deg]','yw [deg]'},'location','EastOutside')
set(gca,'position',axespos)
title('Motion parameters')

ylabel('mm | deg')

ax1=subplot(6,7,7+[1:6]);
x1=[1:size(rp,1)]; y1=rp(:,7); y2=rp(:,8)/pi*180;
[AX,~,~]=plotyy(x1,y1, x1,y2);
set(AX,'xlim',[1,size(rp,1)]);
ylabel(AX(1),'||dShift/dt||_2 (mm)');
ylabel(AX(2),'||dRot/dt||_2 (deg)');
title('L2norm of first derivatives of motion parameters')
xlabel('TR')

sumstat.meanz=zeros(5,size(y,1));
sumstat.meanabsz=zeros(5,size(y,1));
fy = myy_filter(y, TR_sec, FilterBand);
zyres = zscore(fy)';

subplot(6,7,7*2+[1:6]); 
imagesc(zyres);
sumstat.meanz(1,:) = mean(zyres,1);
sumstat.meanabsz(1,:) = mean(abs(zyres),1);
caxis([-2 2]);
hc=colorbaraxes('location','EastOutside');
caxis([-2 2]);
title(hc,'z-score');

R={};
V={};
R{1}=corr(fy);
V{1}=triuval(R{1});
V{1}(isnan(V{1}))=eps;

% setting user-defined regressors to compare
M=[];
if ~isfield(JOB,'covset')
 JOB.covset=[1 2 3 4];
 if isfield(JOB,'nogs')
  JOB.covset=[1 2 4 3];
 end
 % or [1 2 4 3] for rp+cc+art+gs
end
cov_vals={[ones(d(4),1), linspace(-1,1,d(4))',rp], cc, gs};
cov_names={'trend+rigidmotion','+compcor','+globalsignal'};
cov_names_short={'+td+rp',['+cc_n',num2str(JOB.num_pcs)],...
 ['+gs']};
cov_names_shorter={'+td+rp','+cc','+gs'};
tag1{1}='orig';
tag2{1}='orig';
for j=1:3
 COV{j}  = cov_vals{JOB.covset(j)};
 Mdesc{j} = cov_names{JOB.covset(j)};
 tag1{j+1} = cov_names_short{JOB.covset(j)};
 tag2{j+1} = cov_names_shorter{JOB.covset(j)};
end
output_suffix_short = output_suffix;
output_suffix=[output_suffix [tag2{:}]];
title(ax0,[JOB.dir_data,'/',JOB.name_epi,':',output_suffix(2:end)],'interp','none');
sumstat.tag=tag1;

for i=1:3
 M=double([M COV{i}]);
 yres = myy_filter(y-M*(pinv(M'*M)*M'*y), TR_sec, FilterBand);
 subplot(6,7,7*(i+2)+[1:6]);
 zyres=zscore(yres)';
 imagesc(zyres);
 sumstat.meanz(1+i,:) = mean(zyres,1);
 sumstat.meanabsz(1+i,:) = mean(abs(zyres),1);
 caxis([-2 2]); %hc=colorbar;
 title(hc,'z-score');
 if i==3, xlabel('TR'); end
 if i==1, ylabel(['Subsampled voxels with GM>0.95']); end
 
 subplot(6,7,7*(i+3));
 imagesc(zscore(M));
 set(gca,'fontsize',10)
 title(Mdesc{i})
 R{i+1}=corr(yres+eps);
 V{1+i}=triuval(R{i+1});
 V{1+i}(isnan(V{1+i}))=eps;
end

colormap(sgcolormap('CKM'));
name_figure='_timeseries.png';
export_fig([JOB.dir_func,'/',name1,output_suffix,name_figure],['-r',num2str(fig_dpi)])
close(hf);
if isfield(JOB,'plotuntil')&&JOB.plotuntil==1
 return
end

%% plot#2: distribution
isPairDistance=1;
if ~strcmp(version,'7.9.0.529 (R2009b)')
 try 
   pdist
 catch ME
  isPairDistance=0;
 end
end
if isPairDistance
 hf=figure('position',[1 33 933 1061],'visible','off');
 subplot(5,3,[1,2]);
 estker=zeros(4,256);
 Momenta=zeros(4,4);
 for i=1:4
  xi=linspace(-1,1,256);
  estker(i,:) = ksdensity(V{i}, xi);
  Momenta(i,1) = mean(V{i});
  Momenta(i,2) = var(V{i});
  Momenta(i,3) = skewness(V{i});
  Momenta(i,4) = kurtosis(V{i});
 end
 sumstat.momenta=Momenta;
 hold on;
 ColorOrder=[1 0.5 0.5; 0.5 1 0.5; 0.5 0.5 1; 1 0.5 1; 0.5 1 1];
 colororder='rgbmc';
 for i=1:4
  h_lines(i)=plot(xi,estker(i,:),'linewidth',7-i,'color',colororder(i));
 end
 tag=tag1;
 h=legend(h_lines,tag, 'location','eastoutside');
 grid on; box on; set(h,'interp','none')
 title([JOB.dir_data,'/',JOB.name_epi,':',output_suffix_short(2:end)],'interp','none');
 xlabel('Pearson correlation coefficient');ylabel({'estimated','kernel density'})
 
 tag=tag2;
 panelk=[4 7 10 13];
 MomentumDesc={'mean','variance','skewness','kurtosis'};
 for k=2:4
  subplot(5,3,panelk(k)); hold on;
  for j=1:4
   h=barh(j,Momenta(j,k),'facecolor',ColorOrder(j,:));
  end
  set(gca,'ytick',1:5, 'yticklabel',tag, 'ydir','rev');
  xlabel(MomentumDesc{k})
  box on; grid on;
 end
 subplot(5,3,panelk(1));
 barh([max(rp)],'facecolor',[.7 .7 .7]);
 yticklabel={'dx','dy','dz','rx','ry','rz','|dmdt|'};
 title('max')
 set(gca,'ydir','rev','yticklabel',yticklabel)
 xlabel('mm'); grid on; box on;
 
 % correlation over distance
 [ijk1,ijk2,ijk3]=ind2sub(size(gm.img),idx);
 ijk=[ijk1,ijk2,ijk3];
 xyz=ijk2xyz(ijk,gm);
 D=pdist(xyz,'euclidean');
 hold on;
 paneli=[5 8 11 14];
 cmapnames={'WG','WB','WM','WC'};
 for i=1:3
  subplot(5,3,paneli(i))
  dr = V{i+1}-V{i};
  jhist(D,dr,[80 80]);
  [B,~,~,~,stats] = regress(dr,[D', D'*0+1]);
  sumstat.slope(i+1)=B(1);
  hold on;
  plot([0:150], B(2)+B(1)*[0:150],'color','w')
  ylim0=[-max(abs(ylim)), max(abs(ylim))];
  text(5,ylim0(2)*0.8,['b= ',num2str(B(1),2),', p=',num2str(stats(3))]);
  text(5,ylim0(1)*0.8,['R^2= ',num2str(stats(1))]);
  xlim([0 max(D)]); ylim(ylim0);
  xlabel('Euclidean distance(mm)');
  ylabel({'change of corr',['[',tag2{i+1},'] vs. [',tag2{i},']']});
  grid on; box on;
  colormap(sgcolormap(cmapnames{i}));
  freezeColors;
 end
 tag=tag2;
 for i=1:4
  subplot(5,3,i*3)
  h=qqplot(V{i});
  title('');
  set(h(1),'MarkerEdgeColor',ColorOrder(i,:))
  set(h(3),'color','k','linestyle','-');
  [~,p]=kstest(zscore(V{i}));
  sumstat.kspval(i)=p;
  text(-4.7,0.8,['ks:gof=',num2str(p*100,3),'%']);
  grid on; box on;xlabel('norm-q');
  ylabel([tag{i},'-q'],'interp','none');
  axis([-5 5 -1 1]);
 end
 
 % save summary file
 fname_sum=[path1,name1,output_suffix,'.mat'];
 save(fname_sum, 'sumstat');
 
 % save figure
 name_figure='_corrdist.png';
 screen2png([JOB.dir_func,'/',name1,output_suffix,name_figure],fig_dpi);
 close(hf);
end

if isfield(JOB,'plotuntil')&&JOB.plotuntil==2
 return
end

%% plot#3: DMN slice maps
sc=2.5; % figure scale

if isfield(JOB,'dir_fs')
 hf=figure('position',[1972, 33, d(2)*4*sc, d(1)*5*sc], 'color','k');
 %t1w = load_uns_nii([JOB.dir_data,'/oBrain.nii']);
 t1w = load_uns_nii([JOB.dir_data,'/obmt1w.nii']);
 if JOB.sag
  % get a sagittal slice
  t1w = squeeze(t1w.img(crrmap.ijk(1),:,:))';
  crrmap.ind_gm=find(gm.img(crrmap.ijk_med(1),:,:));
 else
  % get an axial slice
  t1w = t1w.img(:,:,crrmap.ijk(3))';
  crrmap.ind_gm=find(gm.img(:,:,crrmap.ijk_med(3)));
 end
 [~,ax] =getLayout(20,[5,4]);
 [~,ax2] =getLayout(10,[5,2],[0.15 0.15]);
 
 M=[];
 for i=1:5
  if i == 1
   crrmap.yres = myy_filter(crrmap.y, TR_sec, FilterBand);
  else
   M=double([M COV{i-1}]);
   %crrmap.yres = myy_filter(crrmap.y-M*((M'*M)\M'*crrmap.y), TR_sec, FilterBand);
   crrmap.yres = myy_filter(crrmap.y-M*(pinv(M'*M)*M'*crrmap.y), TR_sec, FilterBand);
  end
  
  % spatial smoothing
  if JOB.sag
   % 3D (x,y,t)
   crrmap.yres_3d = reshape(crrmap.yres',d(2),d(3),[]);
  else
   crrmap.yres_3d = reshape(crrmap.yres',d(1),d(2),[]);
  end
  numScans=size(crrmap.yres_3d,3);
  for t=1:numScans
   crrmap.yres_3d(:,:,t) = gaussblur(crrmap.yres_3d(:,:,t),1.5);
  end
  crrmap.yres = reshape(crrmap.yres_3d,[],numScans)';
  
  % correlation with mean of precuneous mask on the slice
  crrmap.rsfc = corr(crrmap.yres, mean(crrmap.yres(:,crrmap.idx_prec),2));
  
  % axial slices
  h=axespos(ax,(i-1)*4+4);
  img1=double(t1w);
  if JOB.sag
   img2=reshape(crrmap.rsfc', [d(2) d(3)])';
  else
   img2=reshape(crrmap.rsfc', [d(1) d(2)])';
  end
  imagecorr1(img1, img2, 0.25);
  hold on;
  if JOB.sag
   line([crrmap.ijk_med(2);crrmap.ijk_med(2)],[0 d(3)],  'color','w');
   line([0 d(2)], [crrmap.ijk_med(3);crrmap.ijk_med(3)], 'color','w');
   if i==1, text(d(3)*0.3,d(2)*0.95,'FWHM=1.5 pixels', 'color','w', 'fontsize',12); end
  else
   line([crrmap.ijk_med(1);crrmap.ijk_med(1)],[0 d(2)],  'color','w');
   line([0 d(1)], [crrmap.ijk_med(2);crrmap.ijk_med(2)], 'color','w');
   if i==1, text(d(2)*0.3,d(1)*0.95,'FWHM=1.5 pixels', 'color','w', 'fontsize',12); end
   text(d(1)*0.07, d(2)*0.07, 'R','color','w','fontsize',12)
  end
  
  % and z-scored timeseries
  h=axespos(ax2,(i-1)*2+1);
  if JOB.sag
   imagesc( zscore(crrmap.yres(:,crrmap.ind_gm))' );
  else
   imagesc( zscore(crrmap.yres(:,crrmap.ind_gm))' );
  end
  caxis([-2 2]); hc=colorbar; title(hc,'z-score','color','w');
  set(gca,'xcolor','w','ycolor','w')
  tag=tag1;
  title(tag{i},'color','w','fontsize',14,'interp','none')
  if i==3, ylabel(['GM>0.95 voxles in the slice'], ...
    'color','w','fontsize',14); end
  if i==5
   set(gca,'xtick',[]);
   hxlabel=xlabel('TR','color','w','fontsize',14);
   dd=size(crrmap.yres(:,crrmap.ind_gm));
   set(hxlabel,'position',[dd(1)/2, dd(2)*1.05, 1])
  end
  
  % correlation matrix
  h=axespos(ax,(i-1)*4+3);
  imagesc(corr(crrmap.yres(:,crrmap.ind_gm))); axis image;
  caxis([-.5 .5]); hc=colorbar; title(hc,'corr','color','w');
  set(gca,'xcolor','w','ycolor','w')
  if i==5
   set(gca,'xtick',[]);
   hxlabel=xlabel(['GM>0.95 voxels'],'color','w','fontsize',14);
   dd=repmat(size(crrmap.yres(:,crrmap.ind_gm),2),[1,2]);
   set(hxlabel,'position',[dd(1)/2, dd(2)*1.05, 1])
  end
 end
 
 h=axes('position',[0.75 0 0.25 0.02]);
 hb=colorbar('peer',h,'location','South');
 colormap(sgcolormap('BLUE-RED'));
 axis off; caxis([1,7]);
 set(hb,'xtick',1:7,'xticklabel',{'-1','-0.75','-0.5','|0.25|','0.5','0.75','1'});
 cbfreeze(hb);
 
 colormap(sgcolormap('CKM'));
 name_figure='_corrmap.png';
 screen2png([JOB.dir_func,'/',name1,output_suffix,name_figure], fig_dpi);
 close(hf);
end
if isfield(JOB,'plotuntil')&&JOB.plotuntil==3
 return
end

%% plot#4. sample timeseries (combine with plot#3?)
if isfield(JOB,'dir_fs')
 % 1. find voxel indices for precuneous, med-prefrontal, precentral,
 aparc = load_uns_nii(fname_aparc);
 %aseglabels=[12130 11130]; % precuneous
 %aseglabels=[12109 11109]; % post cingulate gyrus
 % 11146 left central sulcus
 % 12146 right central sulcus
 % 11107 12107 mid-anterior cingulate gyrus & sulcus
 % 11116 12116 superior frontal gyrus
 roiname={'Precuneous','Central sulcus','Anterior cingulate gyrus'};
 AsegLabels={[11130, 12130], [11146], [11107, 12107]}; %, [11116, 12116]};
 %ind1 = crrmap.ind; % seed voxel
 roi=[];
 roi(1).ijk = crrmap.ijk_med;
 roi(1).ind = crrmap.ind;
 l=2;
 for li=2:numel(AsegLabels)
  mask = ismember(aparc.img, AsegLabels{li});
  % z-coordinate should be on the same plane as z = crrmap.ijk_med(3)
  mask_sub = find3(mask);
  if ~numel(mask_sub(mask_sub(:,3)==crrmap.ijk_med(3)))
   warning('no voxel on this slice!');
   continue;
  else
   mask_sub = mask_sub(mask_sub(:,3)==crrmap.ijk_med(3),:);
   roi(l).ijk = round(median(mask_sub));
   if ~ismember(aparc.img(roi(l).ijk(1),roi(l).ijk(2),roi(l).ijk(3)), AsegLabels{l})
    % if the median is not in the label, find nearest precuneous voxel
    roi(l).dist = l2norm(mask_sub - repmat(roi(l).ijk,[size(mask_sub,1),1]));
    [~,ii] = min(roi(l).dist);
    roi(l).ijk = mask_sub(ii,:);
    % double check
    if ~ismember(aparc.img(roi(l).ijk(1),roi(l).ijk(2),roi(l).ijk(3)), AsegLabels{l})
     error('??!!')
    end
   end
   roi(l).ind = sub2ind(d(1:2), roi(l).ijk(1), roi(l).ijk(2)); % index for 2-D slice
   l=l+1;
  end
 end
 
 hf=figure('position',[1972, 33, d(2)*4*sc, d(1)*5*sc], 'color','k');
 t1w = load_uns_nii([JOB.dir_data,'/oBrain.nii']);
 t1w = t1w.img(:,:,crrmap.ijk(3))';
 [~,ax] =getLayout(20,[5,4]);
 [~,ax2] =getLayout(20,[5,4],[0.10 0.15]);
 crrmap.ind_gm=find(gm.img(:,:,crrmap.ijk_med(3)));
 
 M=[];
 for i=1:5
  if i == 1
   crrmap.yres = myy_filter(crrmap.y, TR_sec, FilterBand);
  else
   M=double([M COV{i-1}]);
   %crrmap.yres = myy_filter(crrmap.y-M*((M'*M)\M'*crrmap.y), TR_sec, FilterBand);
   crrmap.yres = myy_filter(crrmap.y-M*(pinv(M'*M)*M'*crrmap.y), TR_sec, FilterBand);
  end
  
  % spatial smoothing!
  crrmap.yres_3d = reshape(crrmap.yres',d(1),d(2),[]);
  numScans=size(crrmap.yres_3d,3);
  for t=1:numScans
   crrmap.yres_3d(:,:,t) = gaussblur(crrmap.yres_3d(:,:,t),1.5);
  end
  crrmap.yres = reshape(crrmap.yres_3d,[],numScans)';
  % correlation with mean of precuneous mask on the slice
  crrmap.rsfc = corr(crrmap.yres, mean(crrmap.yres(:,crrmap.idx_prec),2));
  
  % axial slice: correlation map
  h=axespos(ax,(i-1)*4+1);
  img1=double(t1w);
  img2=reshape(crrmap.rsfc', [d(1) d(2)])';
  if ~exist('corrthres','var')
   corrthres=0.25;
  end
  imagecorr1(img1, img2, corrthres);
  hold on;
  line([crrmap.ijk_med(1);crrmap.ijk_med(1)],[0 d(2)], 'color','w');
  line([0 d(1)], [crrmap.ijk_med(2);crrmap.ijk_med(2)], 'color','w');
  if i==1, text(d(2)*0.3,d(1)*0.95,'FWHM=1.5 pixels', 'color','w', 'fontsize',12); end
  colororder='gmy';
  for l=1:numel(roi)
   scatter(roi(l).ijk(1), roi(l).ijk(2), 100, ...
    'MarkerFaceColor', colororder(l), 'lineWidth',2, 'MarkerEdgeColor','k');
  end
  text(d(1)*0.07, d(2)*0.07, 'R','color','w','fontsize',12)
  
  % individual timeseries
  xlim1=1;      xlim2=min(150, size(crrmap.yres,1));
  if isfield(JOB,'xlim')
   xlim1=JOB.xlim(1);    xlim2=JOB.xlim(2);
  end
  for l=1:numel(roi)
   h=axespos(ax2,(i-1)*4+l+1);
   plot(crrmap.yres(:, roi(l).ind),'color',colororder(l))
   xlim([xlim1 xlim2])
   set(gca,'color','k','xcolor','w','ycolor','w');
   hold on;
   ylim1=ylim;
   text(xlim1+5,ylim1(2)-diff(ylim1)*0.07, roiname{l}, 'color',colororder(l))
   if l==1
    title(tag{i},'fontsize',15,'color','w','interp','none');
   end
   grid on;
   if i==1 && l==1,  xlabel('TR','fontsize',9); end
  end
 end
 
 h=axes('position',[0 0 0.25 0.02]);
 hb=colorbar('peer',h,'location','South');
 colormap(sgcolormap('BLUE-RED'));
 axis off; caxis([1,7]);
 set(hb,'xtick',1:7,'xticklabel',{'-1','-0.75','-0.5','|0.25|','0.5','0.75','1'});
 cbfreeze(hb);
 
 colormap(sgcolormap('CKM'));
 name_figure='_corrmap_inditimeseries.png';
 screen2png([JOB.dir_func,'/',name1,output_suffix,name_figure], fig_dpi);
 close(hf);
end

%% final decision
if ~isfield(JOB,'cov_idx')
 warning([JOB.dir_data,': enter JOB.cov_idx to select regressors and save a residual image']);
 return
end
nii = load_uns_nii(JOB.fname_epi);
% datatype=4, bitpix=16
y = single(reshape(nii.img,[],d(4))');
M=[];
for i=1:JOB.cov_idx
 M=single([M COV{i}]);
end
% LSE (y-y_hat)
yres = y-M*(pinv(M'*M)*M'*y) + mean(y); % putting back the offset

% Band-pass filtering
yres = myy_filter(yres, TR_sec, FilterBand);
% datatype: 8=int16(short), 16=single, 64=double
if d(4)>1000
 nii.hdr.dime.bitpix = 8;
 [yres,offset,slope] = single2int16(yres);
 nii.hdr.dime.scl_inter = offset;
 nii.hdr.dime.scl_slope = slope;
else
 yres = yres-min(yres(:)); % to work with SPM, which only takes positive inputs for fMRI data
 nii.hdr.dime.bitpix = 16;
end
nii.hdr.hist.descrip = ['motion artifacts regressed out: ',output_suffix];
nii.img = reshape(yres',d);
out_prefix=[];
pcsnum=num2str(JOB.num_pcs);
if isfield(JOB,'out_prefix')
 out_prefix=JOB.out_prefix;
else
 if JOB.bpf(1)~=0 && isinf(JOB.bpf(2))
  out_prefix=['hpf'];%r',pcsnum];
 elseif JOB.bpf(1)==0 && ~isinf(JOB.bpf(2))
  out_prefix=['lpf'];%r',pcsnum];
 elseif JOB.bpf(1)~=0 && ~isinf(JOB.bpf(2))
  out_prefix=['bpf'];%r',pcsnum];
 end
 if JOB.cov_idx >= 2
  out_prefix=[out_prefix,'r',pcsnum];
 end
end
fname_cov=[path1,out_prefix,JOB.name_epi];
disp(['> saving noise-regressed timeseries in ',fname_cov,' ..']);
save_untouch_nii(nii, fname_cov);
JOB.fname_cov=fname_cov;

if isfield(JOB,'dir_fs')
 setenv('FSLOUTPUTTYPE','NIFTI');
 [~,name1,~] = fileparts_gz (JOB.name_meanepi);
 myunix(['bet ',path1,name1,' ',path1,name1,'_brain -R -S -f 0.3']);
 myunix(['fslmaths ',path1,name1,'_brain -add ',path1,'o',...
  JOB.t1w_suffix,'.nii -bin ',path1,'epimask.nii']);
 unix(['rm -f ',path1,name1,'_brain*']);
 name2=[out_prefix,JOB.name_epi];
 myunix(['fslmaths ',path1,name2,' -mas ',path1,'epimask ',path1,'m',name2]);
end

cd(path0);

end

function [F,ker] = gaussblur(f,FWHM)
%function F = gaussblur(f,FWHM)
%
% f : input 1D series, 2D image or 3D volume
% FWHM: filter size in sampling point (pixel/voxel).
%
% (c) Moo K. Chung, July 2003. March 2004.
% (cc) sgKIM, 2011, 2015
% + now works for 1D, 2D and 3D.
% + kernel size extended to ensure the smoothness as intended by FWHM.
% + relationship between t and FWHM is corrected now.
%
% See Chung et al., 2001. A Unified Statistical Appraoch to
% Deformation-based morphometry, NeuroImage for implementation detail.


% ksize is the number of neighboring voxels we are averaging.
% it depends on FWHM.
% ksize is allowd to be an odd interger.
ksize = round(FWHM*2+1);

% see Chung et al.(2001) for the relationship between t and FWHM.
% FWHM = 4*sqrt(log(2))*sqrt(t)
%
% FWHM = sqrt(8*log(2)) * sigma (mathworld.Wolfram.com)
% sigma = FWHM / sqrt(8*log(2))
% t = (sigma^2)/2
t = (FWHM^2) / (16*log(2));

d = round(ksize/2);

% we truncate and normalize the Gaussian kernel for fast computation
% 1D, 2D, 3D version of kerenl is slightly different
n = ndims(f);
if n==1
 domain_x=(1:ksize)-d;
 ker=exp(-(domain_x.^2)/(4*t)) /sqrt(t/2);  % not in exponential, so doesn't matter due to normalization
elseif n==2
 domain_x=kron(ones(ksize,1),[1:ksize])-d;
 domain_y=kron(ones(1,ksize),[1:ksize]')-d;
 ker=exp(-(domain_x.^2 +domain_y.^2)/(4*t)) /t;
elseif n==3
 domain_x = repmat ( kron(ones(ksize,1),[1:ksize])-d, [1,1,ksize]);
 domain_y = repmat ( kron(ones(1,ksize),[1:ksize]')-d, [1,1,ksize]);
 domain_z = reshape( repmat ( kron(ones(ksize,1),[1:ksize])-d, [ksize, 1, 1]), [ksize,ksize,ksize]);
 ker=exp(-(domain_x.^2 + domain_y.^2 +domain_z.^2)/(4*t)) /t^(3/2);
end;
ker=ker/sum(ker(:));

% smoothing is done by convolution
F=convn(f,ker,'same');

end
