function JOB = myspm_denoise_slab(JOB)
% JOB = myspm_denoise_slab(JOB)
%
% JOB requires:
%  .name_epi;
%  .cov_idx    1x1 for rp+cc+gs+scrub
%  .num
% (.nofigure)
%
% [example]
% job1=[];
% job1.num_pcs  = 48; (by default; see below for other # of PCs)
% job1.subjid   = 1002;
% job1.dir_fs   = dir1;
% job1.fsd      = 'boldloc';
% job1.name_epi = 'auloc?.nii';
% job1.TR_sec   = 1.5+1.4;
% job1.num_runs = 6;
% job1.dir_base = '/scr/vatikan3/Tonotopy/func/';
% job1.dir_figure = [dir1,'/FIGs/compcor'];
% myspm_compcor_slab(job1)
%
% NOTE: myy_compcor_slab computes 48 PCs by default. For any smaller number of
% PCs, it just takes the first n column vectors.
%
% (cc) 2015-2017, sgKIM.

%% 0. set parameters
path0=pwd;
if ~isfield(JOB,'overwrite'),  overwrite=0; else overwrite=JOB.overwrite; end
subjid = JOB.subjid;
TR_sec = JOB.TR_sec;
if ~isfield(JOB,'num_pcs'), JOB.num_pcs=48;    end;
if ~isfield(JOB,'detrend'), JOB.detrend=1;    end;
if ~isfield(JOB,'varnorm'), JOB.varnorm=1;    end;
if ~isfield(JOB,'bpf1'),    JOB.bpf1=[1/128 Inf]; end;
FilterBand = JOB.bpf1;
if ~isfield(JOB,'param_cc')
JOB.param_cc = sprintf('_n%df%0.2f-%0.2f', JOB.num_pcs, JOB.bpf1);
end
sumstat=[];
if isfield(JOB,'gm_thrs_std') ...
&& isfield(JOB,'mov_thrs_mm') && isfield(JOB,'mov_thrs_deg')
if ~isfield(JOB,'gm_thrs_std'),  JOB.gm_thrs_std  = 5;    end
if ~isfield(JOB,'mov_thrs_mm'),  JOB.mov_thrs_mm  = 0.5;    end
if ~isfield(JOB,'mov_thrs_deg'), JOB.mov_thrs_deg = 0.005; end
mac_output_suffix=sprintf('_%0.1fstd_%0.2fmm_%0.2fdeg', ...
JOB.gm_thrs_std, JOB.mov_thrs_mm, JOB.mov_thrs_deg);
else
mac_output_suffix='_5.0std_0.50mm_0.01deg';
end
if ~isfield(JOB,'t1w_suffix')
JOB.t1w_suffix='t1w';
end
idx=strfind(JOB.name_epi,'?');
switch JOB.name_epi(1:2)
case {'au', 'ua'}
fnameprefix_func = JOB.name_epi(3:idx-1);
case {'a', 'u'}
fnameprefix_func = JOB.name_epi(2:idx-1);
otherwise
fnameprefix_func = JOB.name_epi(1:idx-1);
end

% directory for figures
if ~isfield(JOB,'dir_figure')
dir_figure = [JOB.dir_base,'/',subjid,'/fig_denoise'];
else
dir_figure = JOB.dir_figure;
end
JOB.dir_figure = dir_figure;
[~,~]=mkdir(dir_figure);

% gray matter masks
for r=1:JOB.num_runs
JOB.fname_gmmask{r}=myls([JOB.dir_base,'/',subjid,'/fsc1UNI_run',num2str(r),'.nii']);
end

%% I. Compcor regressors
fname_eigvec = ['cc_',fnameprefix_func,num2str(JOB.num_runs), ...
'_eigvec',JOB.param_cc,'.txt'];
if ~isfield(JOB,'overwrite_cc'), JOB.overwrite_cc=[0 0 0]; end
if ~exist(fname_eigvec,'file') || sum(JOB.overwrite_cc)
job1 = JOB;
job1.subjid =subjid;
job1.overwrite = JOB.overwrite_cc;
myy_compcor_slab(job1); % creating cc_*_gs.txt, cc_*_eigval*.txt cc_*_eigvec*.txt
end

%% II.  MAC: find outliers
for r=1:JOB.num_runs
idx=strfind(JOB.name_epi,'?');
name_epi = [JOB.name_epi(1:idx-1),num2str(r),JOB.name_epi(idx+1:end)];
path1=[fullfile(JOB.dir_base,subjid),'/'];
JOB.fname_epi = fullfile(path1, name_epi);
%path2=[fullfile(JOB.dir_base,subjid),'/fig_denoise/'];
cd(path1);
[~,~]=mkdir(dir_figure);
nii = load_uns_nii(JOB.fname_epi);
d = nii.hdr.dime.dim(2:5);
y = double(img2y(nii.img));
gm = load_uns_nii(JOB.fname_gmmask{r});
y = y(:,gm.img(:)>0.99);

if JOB.num_pcs >0
cc = load([path1,fname_eigvec]);
else
cc = [];
end
% why is it "mac"?
fname_mov = [path1,'/mac_mov_',fnameprefix_func,num2str(r), ...
mac_output_suffix,'.txt']; % rigid motion parameter with length of derivatives
fname_out = [path1,'/mac_out_',fnameprefix_func,num2str(r), ...
mac_output_suffix,'.txt']; % outling frames (scrubbing regressors)

if ~exist(fname_mov,'file') || ~exist(fname_out,'file') || overwrite
job1=JOB;
job1.fname_rp = [path1,'rp_',fnameprefix_func,num2str(r),'.txt'];
myspm_mac(job1);
end

rp = load(fname_mov);
rmparam = load(fname_out);
fname_gs = [path1,'cc_gm99.txt'];
if exist('fname_gs','file') || overwrite
gs = load(fname_gs);
else
gs = ones(size(rmparam,1),1);
end
output_suffix = [JOB.param_cc,'_',mac_output_suffix(2:end), ...
'_b',num2str(FilterBand(1),2),'-',num2str(FilterBand(2),2)];

%% III. FIGURES
if ~isfield(JOB,'nofigures') && (JOB.num_pcs<24)
%% plot#1: timeseries
hf=figure('position',[1952        -172         962        1090]);
ax1=subplot(6,7,[1:6]); hold on;
plot(rp(:,end-1),'r');
plot(rp(:,end)*100,'g');
xlim([1,d(4)]); h=colorbar; set(h,'visible','off');
ylabel({'Frame-wise','disp. [mm|100*deg]'});

sumstat.meanz=zeros(5,size(y,1));
sumstat.meanabsz=zeros(5,size(y,1));
fy = myy_filter(y, TR_sec, FilterBand);
zyres = zscore(fy)';

subplot(6,7,7+[1:6]); imagesc(zyres);
sumstat.meanz(1,:) = mean(zyres,1);
sumstat.meanabsz(1,:) = mean(abs(zyres),1);
caxis([-2 2]); hc=colorbar; title(hc,'z-score');

R={};
V={};
rndidx = randi(size(fy,2),[256,1]);
R{1}=corr(fy(:,rndidx));
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
cov_vals={[ones(d(4),1), linspace(-1,1,d(4))',rp], cc, gs, rmparam};
cov_names={'trend+rigidmotion','+compcor','+globalsignal','+scrubbing'};
cov_names_short={'+td+rp',['+cc_th0.99_n',num2str(JOB.num_pcs)],...
['+gs_0.99'],['+scrb_',num2str(size(rmparam,2))]};
cov_names_shorter={'+td+rp','+cc','+gs','+scrb'};
tag1{1}='orig';
tag2{1}='orig';
for j=1:4
COV{j}  = cov_vals{JOB.covset(j)};
Mdesc{j} = cov_names{JOB.covset(j)};
tag1{j+1} = cov_names_short{JOB.covset(j)};
tag2{j+1} = cov_names_shorter{JOB.covset(j)};
end
output_suffix_short = output_suffix;
output_suffix=[output_suffix [tag2{:}]];

title(ax1,[subjid,'/',JOB.name_epi,':',output_suffix(2:end)],'interp','none');
sumstat.tag=tag1;

for i=1:4
M=double([M COV{i}]);
yres = myy_filter(y-M*(pinv(M'*M)*M'*y), TR_sec, FilterBand);
subplot(6,7,7+7*i+[1:6]);
zyres=zscore(yres)';
imagesc(zyres);
sumstat.meanz(1+i,:) = mean(zyres,1);
sumstat.meanabsz(1+i,:) = mean(abs(zyres),1);
caxis([-2 2]); hc=colorbar; title(hc,'z-score');
if i==4, xlabel('TR'); end
if i==2, ylabel(['Subsampled voxels with GM>0.99']); end

subplot(6,7,7+7*i+7);
imagesc(zscore(M));
set(gca,'fontsize',10)
title(Mdesc{i})
R{i+1}=corr(yres(:,rndidx)+eps);
V{1+i}=triuval(R{i+1});
V{1+i}(isnan(V{1+i}))=eps;
end

colormap(sgcolormap('CKM'));
name_figure='_timeseries.png';
screen2png([dir_figure,'/resy',output_suffix,name_figure]);
close(hf);

end % if-figure
end % run-loop
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


%% old codes for visualization
%      if isfield(JOB,'plotuntil')&&JOB.plotuntil==1
%       continue
%      end
%
%    %% plot#2: distribution
%    hf=figure('position',[1972 33 933 1061]);
%    subplot(5,3,[1,2]);
%    estker=zeros(5,256);
%    Momenta=zeros(5,4);
%    for i=1:5
%     xi=linspace(-1,1,256);
%     estker(i,:) =ksdensity(V{i}, xi);
%     Momenta(i,1) = mean(V{i});
%     Momenta(i,2) = var(V{i});
%     Momenta(i,3) = skewness(V{i});
%     Momenta(i,4) = kurtosis(V{i});
%    end
%    sumstat.momenta=Momenta;
%    hold on;
%    ColorOrder=[1 0.5 0.5; 0.5 1 0.5; 0.5 0.5 1; 1 0.5 1; 0.5 1 1];
%    colororder='rgbmc';
%    for i=1:5
%     h_lines(i)=plot(xi,estker(i,:),'linewidth',7-i,'color',colororder(i));
%    end
%    tag=tag1;
%    h=legend(h_lines,tag, 'location','eastoutside');
%    grid on; box on; set(h,'interp','none')
%    title([subjid,'/',JOB.name_epi,':',output_suffix_short(2:end)],'interp','none');
%    xlabel('Pearson correlation coefficient');ylabel({'estimated','kernel density'})
%
%    tag=tag2;
%    panelk=[4 7 10 13];
%    MomentumDesc={'mean','variance','skewness','kurtosis'};
%    for k=2:4
%     subplot(5,3,panelk(k)); hold on;
%     for j=1:5
%      h=barh(j,Momenta(j,k),'facecolor',ColorOrder(j,:));
%     end
%     set(gca,'ytick',1:5, 'yticklabel',tag, 'ydir','rev');
%     xlabel(MomentumDesc{k})
%     box on; grid on;
%    end
%    subplot(5,3,panelk(1));
%    barh([max(rp)],'facecolor',[.7 .7 .7]);
%    yticklabel={'dx','dy','dz','rx','ry','rz','FDt','FDr'};
%    title('max')
%    set(gca,'ydir','rev','yticklabel',yticklabel)
%    xlabel('mm'); grid on; box on;
%
%    tag=tag2;
%    for i=1:5
%     subplot(5,3,i*3)
%     h=qqplot(V{i});
%     title('');
%     set(h(1),'MarkerEdgeColor',ColorOrder(i,:))
%     set(h(3),'color','k','linestyle','-');
%     [~,p]=kstest(zscore(V{i}));
%     sumstat.kspval(i)=p;
%     text(-4.7,0.8,['ks:gof=',num2str(p*100,3),'%']);
%     grid on; box on;xlabel('norm-q');
%     ylabel([tag{i},'-q'],'interp','none');
%     axis([-5 5 -1 1]);
%    end
%
%    % save summary file
%    fname_sum=[path1,'resy',output_suffix,'.mat'];
%    save(fname_sum, 'sumstat');
%
%    % save figure
%    name_figure='_corrdist.png';
%    screen2png([path2,'resy',output_suffix,name_figure]);
%    close(hf);
%    if isfield(JOB,'dir_figure')
%     copyfile([path2,'resy',output_suffix,name_figure],...
%      [JOB.dir_figure,'/resy',output_suffix,name_figure(1:end-4),'_',subjid,'.png']);
%    end
%
%    if isfield(JOB,'plotuntil')&&JOB.plotuntil==2
%     continue
%    end
%    %% plot#4. sample timeseries (combine with plot#3?)
%    if isfield(JOB,'sampletimeseries');
%     % 1. find voxel indices for precuneous, med-prefrontal, precentral,
%     aparc = load_uns_nii(fname_aparc);
%     %aseglabels=[12130 11130]; % precuneous
%     %aseglabels=[12109 11109]; % post cingulate gyrus
%     % 11146 left central sulcus
%     % 12146 right central sulcus
%     % 11107 12107 mid-anterior cingulate gyrus & sulcus
%     % 11116 12116 superior frontal gyrus
%     roiname={'Precuneous','Central sulcus','Anterior cingulate gyrus'};
%     AsegLabels={[11130, 12130], [11146], [11107, 12107]}; %, [11116, 12116]};
%     %ind1 = crrmap.ind; % seed voxel
%     roi=[];
%     roi(1).ijk = crrmap.ijk_med;
%     roi(1).ind = crrmap.ind;
%     l=2;
%     for li=2:numel(AsegLabels)
%      mask = ismember(aparc.img, AsegLabels{li});
%      % z-coordinate should be on the same plane as z = crrmap.ijk_med(3)
%      mask_sub = find3(mask);
%      if ~numel(mask_sub(mask_sub(:,3)==crrmap.ijk_med(3)))
%       warning('no voxel on this slice!');
%       continue;
%      else
%       mask_sub = mask_sub(mask_sub(:,3)==crrmap.ijk_med(3),:);
%       roi(l).ijk = round(median(mask_sub));
%       if ~ismember(aparc.img(roi(l).ijk(1),roi(l).ijk(2),roi(l).ijk(3)), AsegLabels{l})
%        % if the median is not in the label, find nearest precuneous voxel
%        roi(l).dist = l2norm(mask_sub - repmat(roi(l).ijk,[size(mask_sub,1),1]));
%        [~,ii] = min(roi(l).dist);
%        roi(l).ijk = mask_sub(ii,:);
%        % double check
%        if ~ismember(aparc.img(roi(l).ijk(1),roi(l).ijk(2),roi(l).ijk(3)), AsegLabels{l})
%         error('??!!')
%        end
%       end
%       roi(l).ind = sub2ind(d(1:2), roi(l).ijk(1), roi(l).ijk(2)); % index for 2-D slice
%       l=l+1;
%      end
%     end
%
%     sc=2.5;
%     corrthres=0.25;
%     hf=figure('position',[1972, 33, d(2)*4*sc, d(1)*5*sc], 'color','k');
%     t1w = load_uns_nii([JOB.dir_base,subjid,'/oBrain.nii']);
%     t1w = t1w.img(:,:,crrmap.ijk(3))';
%     [~,ax] =getLayout(20,[5,4]);
%     [~,ax2] =getLayout(20,[5,4],[0.10 0.15]);
%     crrmap.ind_gm=find(gm.img(:,:,crrmap.ijk_med(3)));
%
%     M=[];
%     for i=1:5
%      if i == 1
%       crrmap.yres = myy_filter(crrmap.y, TR_sec, FilterBand);
%      else
%       M=double([M COV{i-1}]);
%       crrmap.yres = myy_filter(crrmap.y-M*((M'*M)\M'*crrmap.y), TR_sec, FilterBand);
%      end
%
%      % spatial smoothing!
%      crrmap.yres_3d = reshape(crrmap.yres',d(1),d(2),[]);
%      numScans=size(crrmap.yres_3d,3);
%      for t=1:numScans
%       crrmap.yres_3d(:,:,t) = gaussblur(crrmap.yres_3d(:,:,t),1.5);
%      end
%      crrmap.yres = reshape(crrmap.yres_3d,[],numScans)';
%      % correlation with mean of precuneous mask on the slice
%      crrmap.rsfc = corr(crrmap.yres, mean(crrmap.yres(:,crrmap.idx_prec),2));
%
%      % axial slice: correlation map
%      h=axespos(ax,(i-1)*4+1);
%      img1=double(t1w);
%      img2=reshape(crrmap.rsfc', [d(1) d(2)])';
%      imagecorr1(img1, img2, corrthres);
%      hold on;
%      line([crrmap.ijk_med(1);crrmap.ijk_med(1)],[0 d(2)], 'color','w');
%      line([0 d(1)], [crrmap.ijk_med(2);crrmap.ijk_med(2)], 'color','w');
%      if i==1, text(d(2)*0.3,d(1)*0.95,'FWHM=1.5 pixels', 'color','w', 'fontsize',12); end
%      colororder='gmy';
%      for l=1:numel(roi)
%       scatter(roi(l).ijk(1), roi(l).ijk(2), 100, ...
%        'MarkerFaceColor', colororder(l), 'lineWidth',2, 'MarkerEdgeColor','k');
%      end
%      text(d(1)*0.07, d(2)*0.07, 'R','color','w','fontsize',12)
%
%      % individual timeseries
%      xlim1=1;      xlim2=min(150, size(crrmap.yres,1));
%      if isfield(JOB,'xlim')
%       xlim1=JOB.xlim(1);    xlim2=JOB.xlim(2);
%      end
%      for l=1:numel(roi)
%       h=axespos(ax2,(i-1)*4+l+1);
%       plot(crrmap.yres(:, roi(l).ind),'color',colororder(l))
%       xlim([xlim1 xlim2])
%       set(gca,'color','k','xcolor','w','ycolor','w');
%       hold on;
%       ylim1=ylim;
%       text(xlim1+5,ylim1(2)-diff(ylim1)*0.07, roiname{l}, 'color',colororder(l))
%       if l==1
%        title(tag{i},'fontsize',15,'color','w','interp','none');
%       end
%       grid on;
%       if i==1 && l==1,  xlabel('TR','fontsize',9); end
%      end
%     end
%
%     h=axes('position',[0 0 0.25 0.02]);
%     hb=colorbar('peer',h,'location','South');
%     colormap(sgcolormap('BLUE-RED'));
%     axis off; caxis([1,7]);
%     set(hb,'xtick',1:7,'xticklabel',{'-1','-0.75','-0.5','|0.25|','0.5','0.75','1'});
%     cbfreeze(hb);
%
%     colormap(sgcolormap('CKM'));
%     name_figure='_corrmap_inditimeseries.png';
%     screen2png([path2,'resy',output_suffix,name_figure]);
%     close(hf);
%     if isfield(JOB,'dir_figure')
%      copyfile([path2,'resy',output_suffix,name_figure],...
%       [JOB.dir_figure,'/resy',output_suffix,name_figure(1:end-4),'_',subjid,'.png']);
%     end
%    end


%  %% final decision?
%  if ~isfield(JOB,'cov_idx')
%   disp([subjid,': enter JOB.cov_idx to select regressors and save a residual image']);
%   continue
%  else
%   nii = load_uns_nii(fname_epi);
%   % datatype=4, bitpix=16
%   y = single(reshape(nii.img,[],d(4))');
%   M=[];
%   for i=1:JOB.cov_idx
%    M=single([M COV{i}]);
%   end
%   % LSE (y-y_hat) and band-pass filtering
%   yres = myy_filter(y-M*(pinv(M'*M)*M'*y), TR_sec, FilterBand);
%   % don't forget to enter space x time, instead of time x space
%   if  (nii.hdr.dime.bitpix == 16) && (d(4)>1000)
%    % datatype: 8=int16(short), 16=single, 64=double
%    [yres,offset,slope] = single2int16(yres);
%    nii.hdr.hist.descrip = sprintf(['motion artifacts regressed out: %s, ', ...
%     'single2int16:offset=%f, slope=%f'], output_suffix, offset, slope);
%   else
%    nii.hdr.hist.descrip = sprintf(['motion artifacts regressed out: %s'], ...
%     output_suffix);
%   end
%   nii.img = reshape(yres',d);
%   pcsnum=num2str(JOB.num_pcs);
%   if isfield(JOB,'out_prefix')
%    out_prefix=JOB.out_prefix;
%   else
%    out_prefix=['fr',pcsnum];
%   end
%   fname_out=[path1,out_prefix,JOB.name_epi];
%   disp(['> saving residual in ',fname_out,'..']);
%   save_untouch_nii(nii, fname_out);
%  end
