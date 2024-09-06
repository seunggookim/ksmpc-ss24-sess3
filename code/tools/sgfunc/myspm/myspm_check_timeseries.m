function JOB = myspm_check_timeseries (JOB)
% JOB = myspm_check_timeseries (JOB)
%
% JOB
%  .dir_base
%  .fnames_proc
%  .names_proc
%  .fname_art
%  .subjID
%  .dir_png
%
% (cc) 2015, sgKIM.   solleo@gmail.com   https:ggooo.wordpress.com

subjID = fsss_subjID(JOB.subjID);
dir0=pwd;

if ~isfield(JOB,'cthres'), cthres = 0.99; else cthres = JOB.cthres; end

for n=1:numel(subjID)
subjid = subjID{n};
cd(fullfile(JOB.dir_base,subjid));

if ~isfield(JOB,'dir_png')
dir_png = pwd;
else
dir_png = JOB.dir_png;
end
if isfield(JOB,'names_proc')
JOB.fname_png=[dir_png,'/',subjid,'_timecourse_', ...
JOB.names_proc{1},'_vs_',JOB.names_proc{2},'.png'];
else
JOB.fname_png=[dir_png,'/',subjid,'_timecourse_', ...
JOB.fnames_proc{1},'_vs_',JOB.fnames_proc{2},'.png'];  end

C={'gm','wm','cf'};
for c=1:3
[~,fnames] = mydir(['wc',num2str(c),'*.nii']);
nii = load_uns_nii(fnames{1});
IDX{c} = zeroone(nii.img)>cthres;
disp(['# of voxels with Pr(',C{c},')>',num2str(cthres,1),': ', ...
num2str(sum(IDX{c}(:))), ...
' (',num2str(sum(IDX{c}(:))*prod(nii.hdr.dime.pixdim(2:4))/(1000)^2,2),'L)']);
end
clear nii

figure('position',[1921  1   824 976]);

for i=1:2
subplot(5,1,[2*i,2*i+1])
nii = load_uns_nii(JOB.fnames_proc{i});
d = size(nii.img);
index = IDX{1}|IDX{2}|IDX{3};
y=zeros(d(4), sum(index(:)));
for t=1:d(4)
vol = nii.img(:,:,:,t);
y(t,:) = vol(index(:));
end
% tissue class marker
tcm = IDX{1}*0;
tcm(IDX{1}) = -1;
tcm(IDX{2}) = 0;
tcm(IDX{3}) = 1;
[tc,reorder] = sort(tcm(index));
y = y(:,reorder);
%scp = sigchangeperc(y);
%imagesc([repmat(tc,[1,6])'; scp]')
%imagesc([repmat(tc,[1,6])'; zscore(y)]')
imagesc(zscore(y)');
caxis([-5 5]); 
h=colorbar;
%ylabel(h,'Signal change (%)');
ylabel(h,'Z-score','fontsize',12);
colormap(sgcolormap('CKM'));
ylabel(['Voxels with Pr(tissue)>.99'],'fontsize',14);% with Pr(gm)>.99, Pr(wm)>.99, Pr(csf)>.99']);
hold on;
tidx = find(diff(tc));
text(5, 4000, 'GM','color','w','fontsize',14)
line([0 d(4)]',[tidx(1) tidx(1)]-0.5,'color','w','linewidth',2);
text(5, tidx(1)+4000, 'WM','color','w','fontsize',16)
line([0 d(4)]',[tidx(2) tidx(2)]-0.5,'color','w','linewidth',2);
text(5, tidx(2)+4000, 'CSF','color','w','fontsize',16)
if i==2, xlabel('TR','fontsize',12); end;
set(gca,'ydir','nor','xtick',[-5, 0:100:d(4)],'fontsize',12);
if ~isfield(JOB,'names_proc')
title(JOB.fnames_proc{i},'interp','none')
else
title(JOB.names_proc{i},'fontsize',16);
end
end
subplot(5,1,1)
load(JOB.fname_art,'R');
plot(R(:,7)); grid on; box on;
xlim([1 d(4)]);
set(gca,'ydir','nor','fontsize',12);
%xlabel('TR','fontsize',12, 'interp','none');
ylabel(['Mov_Art(mm)'],'fontsize',14, 'interp','none');
h=colorbar; set(h,'visible','off')


screen2png(JOB.fname_png,120);
close(gcf);

end
cd(dir0);
end
