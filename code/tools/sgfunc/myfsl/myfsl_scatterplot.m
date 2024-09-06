function JOB = myfsl_scatterplot (JOB)
%

% JOB requires:
%  .model        [term] SurfStat 'term' format of a GLM
%  .cidx         [1x1] the index of a regressor of interest in the model
% (.numPerm)     [1x1] the number of randomization
%  .dir_base     'Nx1' a directory that encloses a model's subdirectory
%  .meas         'Nx1' measure name 'FA', 'MD', 'L1', "L23", 'GMV', 'WMV', "xmap"
%  .fname_input  'Nx1' fullpath filename of "all_*"
%  .fname_mask   'Nx1' fullpath filename of a mask
% (.glm_prefix)
%
% (cc) 2015, sgKIM.  solleo@gmail.com   https://ggooo.wordpress.com

global overwrite; if isempty(overwrite), overwrite=0; end
if nargin<1, help myfsl_scatterplot; end

%% find the corrp and t-stat files

if ~isfield(JOB,'glm_prefix'), JOB.glm_prefix='';  end
JOB.model_desc = fsss_model_desc(JOB.model, JOB.cidx);
if ~isfield(JOB,'dir_glm')
JOB.dir_glm = [JOB.dir_base,'/GLM/',JOB.meas,'_',JOB.glm_prefix,JOB.model_desc,'/'];
end
M=double(JOB.model);
varnames = char(JOB.model);
if numel(JOB.cidx) == 1
cidx = sort(JOB.cidx);
JOB.vi=varnames{cidx};
x = M(:,cidx);
elseif numel(JOB.cidx) == 2
cidx = JOB.cidx;
JOB.vi = [varnames{cidx(1)},'-',varnames{cidx(2)}];
x = M(:,cidx(1)) - M(:,cidx(2));
end
M(:,cidx) = []; % for covariates
prefix = fullfile(JOB.dir_glm,JOB.vi);
fname_data = fullfile(JOB.dir_glm, ['all_',JOB.meas,'_skeletonised.nii.gz']);
data = load_uns_nii(fname_data);
% thick_pval = load_uns_nii(JOB.fname_thikcpval);
% thick_pval = thick_pval.img;
p1 = load_uns_nii(JOB.fname_pval{1});
p2 = load_uns_nii(JOB.fname_pval{2});
thick_pval = p1.img*0;
thick_pval(~~p1.img) =  p1.img(~~p1.img);
thick_pval(~~p2.img) = -p2.img(~~p2.img);
mean_FA = load_uns_nii(JOB.fname_meanFA);

%% find the local maxima
setenv('FSLOUTPUTTYPE','NIFTI_GZ');
if ~isfield(JOB,'alpha'), JOB.alpha=0.95; end
%sigC=[];
Sign={'+','-'};
for c=1:2
JOB.fname_tfce{c} = [prefix,'_tfce_corrp_tstat',num2str(c),'.nii.gz'];
if myfsl_minp(JOB.fname_tfce{c}) < 0.05
%   sigC=[sigC, c];
JOB.fname_clustxt{c} = [prefix,'_tfce_corrp_tstat',num2str(c),'_clusid.txt'];
if ~exist(JOB.fname_clustxt{c},'file') || overwrite
disp(['# ',JOB.meas,':',JOB.model_desc,'; ',Sign{c},JOB.vi]);
unix(['cluster -i ',JOB.fname_tfce{c},' --thresh=',num2str(JOB.alpha), ...
' --mm > ',JOB.fname_clustxt{c}]);
end
JOB.fname_clusid{c} = [prefix,'_tfce_corrp_tstat',num2str(c),'_clusidmap.nii.gz'];
if ~exist(JOB.fname_clusid{c},'file') || overwrite
unix(['cluster -i ',JOB.fname_tfce{c},' --thresh=',num2str(JOB.alpha), ...
' --mm --oindex=',JOB.fname_clusid{c}]);
end
T = dlmread(JOB.fname_clustxt{c},'\t',1,0);
for j=1:max(T(:,1))
CLUS{c}(j).id = T(j,1);
CLUS{c}(j).numVox = T(j,2);
CLUS{c}(j).pval = 1-T(j,3);

strct = myfsl_atlasquery(T(j,4:6));
CLUS{c}(j).peak_coord_mm = T(j,4:6);
CLUS{c}(j).peak_name = strct.name;
CLUS{c}(j).peak_prob_prcnt = strct.prob;

strct = myfsl_atlasquery(T(j,7:9));
CLUS{c}(j).cog_coord_mm  = T(j,7:9);
CLUS{c}(j).cog_name  = strct.name;
CLUS{c}(j).cog_prob_prcnt = strct.prob;

%% now creat a scatter plot here
ijk = xyz2ijk(CLUS{c}(j).peak_coord_mm, data);
y = squeeze(data.img(ijk(1) , ijk(2), ijk(3), :));
%scatter1(x,y,M, JOB)
hf=figure('position',[ 2103, 296, 640, 306]);
axes('position',[0 0 .5 1]);
ijk = xyz2ijk(CLUS{c}(j).peak_coord_mm, mean_FA);
cfg=struct('ijk', ijk, 'slicedim',3);
if isfield(JOB,'colormap'), cfg.colormap=JOB.colormap; end
baseimg = zeroone(mean_FA.img);
baseimg(baseimg==0) = 1;
imageover1(baseimg, thick_pval, cfg);
set(gca,'xdir','rev')
hold on
x0=xlim; y0=ylim;
line([ijk(1) ijk(1); x0]', [y0; ijk(2) ijk(2)]', 'color','k')
%line([ijk(2) ijk(2); x0]', [y0; ijk(1) ijk(1)]', 'color','k')

%axes('position',[0.5+0.025 0.025 0.45 0.95]);
subplot(122)
[t,df,p]=SurfStatPlot(x,y,M);
set(gca,'xtick',[0 1],'xticklabel',{'NAP','AP'},'fontsize',13);
%vol = num2str(T(j,2)*prod(data.hdr.dime.pixdim(2:4)));
title({CLUS{c}(j).peak_name, sprintf('T(%d)=%0.3f, p=%0.3f, #vox=%d', df,t,p, T(j,2)) },'fontsize',14);
xlim([-0.5 1.5]); ylabel(['Adj. ',JOB.meas],'fontsize',13); grid on;
screen2png(['/scr/vatikan3/APConn/figures/ms_0.52/tbss_',JOB.meas,'_',JOB.model_desc,'_knt',num2str(c),'_clus',pad(j,3),'.png']);
close(hf);
end
end
end
end


%%
function scatter1(x,y,M, JOB)
% find adjusted response variable
[maxT,df,p0,Yadj]=SurfStatPlot(x,y,M, [], 'marker','none', 'markerFaceColor','w', 'LineWidth',2);
%[maxT,~,~,Yadj]=SurfStatPlot(term(x), term(y), term(M), 'marker','none', 'markerFaceColor','w', 'LineWidth',2);
cla; box on; grid on; hold on;
maxT = maxT * (-1+2*(c==1)); % negative for the 2nd contrast
hold on;
Xvalkinds = unique(x);
if numel(Xvalkinds) == 2
% box plots for binary x
h=boxplot(Yadj, x ,'positions', x);
set(findobj(h,'tag','Box'),'color','k');
set(findobj(h,'tag','Median'),'color',[.3 .3 .3]);
set(findobj(h,'tag','Outliers'),'visible','off');
hold on;
else
% regression line for continuous variable
[t,df,p0,Yadj,hs] = SurfStatPlot(x,y,M,[], 'marker','none', ...
'color', [.8 .8 .8], 'LineWidth',2);
h0=hs(1);
set(h0,'linestyle','-');
end
return

% add chance level marker for APS
if strcmpi(JOB.ContrastName{cont_i}(2:end),'APS')
hold on;
% APS = 1 - mace/6
% 1-(mean(abs([0:11]-6)))/6 =0.50
line([.5 .5]', [0 1e+5]','color',[.8 .8 .8],'linewidth',2,'linestyle',':');
axis(CurrentAxis);
end

% using term structure for group labeling
if isfield(JOB,'isGroup')
isG = JOB.isGroup;
if islogical(isG)
idx1 =  ~isG;          idx2 = ~~isG;
else % if it's a term
idx1 =  ~double(isG);  idx2 = ~~double(isG);
end
if isfield(JOB, 'GroupName')
Gtag = JOB.GroupName;
else
if ~islogical(isG)
charG = char(isG);
Gtag = {['Non-',charG{1}], charG{1}};
else
Gtag={'Non-Group','Group'};
end
if JOB.isAsym
Gtag={'','R>L'};
end
end
% open circle for non-Group
sh1=scatter(SLM{kidx,s,cont_i}.X(idx1,knst_i), Yadj(idx1), 8^2, 'o', ...
'markerFaceColor',ones(1,3)*1, 'markerEdgeColor',ones(1,3)*0.25, 'linewidth',1);
% close circle for Group
sh2=scatter(SLM{kidx,s,cont_i}.X(idx2,knst_i), Yadj(idx2), 8^2, 'o','filled', ...
'markerFaceColor',ones(1,3)*.8, 'markerEdgeColor',ones(1,3)*0.25, 'linewidth',1);
% add legend
LegendLocation={'NorthEast','SouthEast','NorthEast'};
lli=1+(sum(JOB.contrasts(cont_i,:))>0)+2*(sum(JOB.contrasts(cont_i,:))<0);
legend([sh1, sh2], Gtag,'location' ,LegendLocation{lli});
else % for no group coding
scatter(SLM{kidx,s,cont_i}.X(:,knst_i), Yadj, 8^2, 'o','filled', ...
'markerFaceColor',ones(1,3)*.8, 'markerEdgeColor',ones(1,3)*0.25, 'linewidth',1);
end

titlestring={txt0, [stattext,'=',sprintf('%0.2f',maxT),...
', area=',sprintf('%0.0f',sum(SURF.AREA{s}(SLM{kidx,s,cont_i}.clusid==clusi))),'mm^2'...
', cluster-P=',sprintf('%0.4f',min(SLM{kidx,s,cont_i}.clus.P(clusi)))]};
title(titlestring,'interp','none','fontsize',14);

if ~isempty(SLM{kidx,s,cont_i}.clusid)
NumLevel=numel(unique(JOB.model(JOB.cidx)));
if NumLevel > 2
xname=JOB.ContrastName{cont_i}(2:end);
if isfield(JOB,'xticklabel') && isfield(JOB,'xtick')  % specifically defined xtick, xticklabel
set(gca,'xtick',JOB.xtick,'xticklabel',JOB.xticklabel);
end
elseif NumLevel == 2
if isfield(JOB,'xticklabel')
set(gca,'xtick',unique(JOB.model(JOB.cidx)),'xticklabel',JOB.xticklabel);
else
xname='Group';
set(gca,'xtick',unique(JOB.model(JOB.cidx)),'xticklabel',{['Non-',xname],xname});
end
end

if numel(unique(SLM{kidx,s,cont_i}.X(:,knst_i))) == 2 % when it's group comparison
xlim([-0.5 1.5]);
%xlabel('Group','fontsize',16);
else
if isempty(x_unit)
xlabel(xname,'fontsize',16);
else
xlabel([xname,' (',x_unit,')'],'fontsize',16);
end
end
if ~isempty(M)
ylabel(['Adj. ',YLabel],'fontsize',16);
else
ylabel(YLabel,'fontsize',16);
end

%% set output jpeg filename
fname_jpeg=fullfile(dir_figure,[JOB.meastype,'_',fname_prefix,'_', ...
JOB.ContrastName{cont_i},'_s',fwhm_mm,'mm',masktxt, ...
'_k',num2str(k),'_s',num2str(s),'_cluster',num2str(clusi),'.jpg']);
if isfield(JOB,'numperms')
[a,b,c] = fileparts(fname_jpeg);
fname_jpeg=[a,'/',b,'_p',num2str(JOB.numperms),c];
end
if isfield(JOB,'fc')
[a,b,c] = fileparts(fname_jpeg);
fname_jpeg=[a,'/',JOB.fc,'_',b,c];
end
disp(['# Saving scatterplots: ',fname_jpeg]);
if isfield(JOB,'figure_dpi')      cfg =[];
screen2jpeg(fname_jpeg, JOB.figure_dpi*2);
else
screen2jpeg(fname_jpeg, 200);
end
close(h_scat);
SLM{kidx,s,cont_i}.H{clusi} = H;
end
end
