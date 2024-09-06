function JOB = myfsl_viewvmap (JOB)
% JOB requires:
%  .dir_base
%  .subjid
%  .dir_sub
%
% (cc) sgKIM 2015

subjID = fsss_subjID(JOB.subjID);
if ~isfield(JOB,'vismethod'),  JOB.vismethod = 'vmap1'; end
if ~isfield(JOB,'bitprefix'),  JOB.bitprefix = 'spm12_1p7'; end
if ~isfield(JOB,'pthres'),     JOB.pthres=0.7; end

for i=1:numel(subjID)
subjid = subjID{i};
v1  = load_uns_nii([JOB.dir_base,'/',subjid,'/',JOB.bitprefix,'_Tensor_V1.nii.gz']);
fa  = load_uns_nii([JOB.dir_base,'/',subjid,'/',JOB.bitprefix,'_Tensor_FA.nii.gz']);
t1w = load_uns_nii([JOB.dir_base,'/',subjid,'/t1w_',JOB.bitprefix,'.nii']);
d = size(t1w.img);
if ~isfield(JOB,'ijks')
JOB.ijks=round(linspace(d(JOB.viewdim)*0.53,d(JOB.viewdim)*0.75,40)); 
end

%% for each seed
seeds = JOB.SEED_NAMES;
for k=1:numel(seeds)
[~,res] = mydir([JOB.dir_base,'/',subjid,'.probtrkX/',JOB.dir_sub, ...
'/vmap_',seeds{k},'*.nii.gz']);
hf=[];
for j=1:numel(res)
vmap = load_uns_nii(res{j});
cfg=[];
if ~isfield(JOB,'alpha'),  JOB.alpha=0.7; end
if ~isfield(JOB,'pthres'), JOB.pthres=300; end
cfg.alpha  = JOB.alpha;
cfg.pthres = JOB.pthres;
[path1,name0,~] = fileparts_gz(res{j});
if JOB.isCoord
% if this is a coordinate read it from the filename
name1 = name0(numel(['vmap_',seeds{k}]):end);
idx = strfind(name1,'_');
cfg.ijk = round([str2double(name1(idx(1)+1:idx(2)-1)), ...
str2double(name1(idx(2)+1:idx(3)-1)), ...
str2double(name1(idx(3)+1:end))]);
end
if cfg.pthres<1
figure;
subplot(211);hist(vmap.img(vmap.img>0),100);
xlabel(['visits']); ylabel(['# of voxels']);
title(res{j},'interp','none');

% Alfred's suggestion
vmap.img = log(vmap.img);
vmap.img(isinf(vmap.img)) = 0;
vmap.img = zeroone(vmap.img);

subplot(212);hist(vmap.img(vmap.img>0),100);
xlabel(['zero-one(log(visits))']); ylabel(['# of voxels']);
title(res{j},'interp','none');
hold on;
ylim0=ylim;
line([cfg.pthres;cfg.pthres],[0;ylim0(2)],'color','r');
screen2jpeg([path1,'/hist_',name0,'.jpg']);
close(gcf);
end
hf(j)=figure;
switch JOB.vismethod
case 'vmap1'
imagevmap1(t1w, v1, fa, vmap, cfg);
name2 = ['orth_',name0];
case 'vmaps'
cfg.ijks = JOB.ijks;
cfg.viewdim = JOB.viewdim;
imagevmaps(t1w, v1, fa, vmap, cfg);
name2 = ['mntg_',name0];
end
set(hf(j),'color','k'); drawnow;
screen2jpeg([path1,'/',name2],90);
if isfield(JOB,'dir_figure')
[~,~] = mkdir(JOB.dir_figure);
copyfile([path1,'/',name2,'.jpg'], [JOB.dir_figure,'/',subjid,'_',name2,'.jpg']);
end
close(hf(j));
end
end
end

end
