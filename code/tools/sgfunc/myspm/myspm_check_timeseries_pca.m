function JOB = myspm_check_timeseries_pca (JOB)
% JOB = myspm_check_timeseries_pca (JOB)
%
% JOB
%  .dir_base
%  .fname_proc
%  .name_proc
%  .subjID
%  .dir_png
%
% (cc) 2015, sgKIM.   solleo@gmail.com   https:ggooo.wordpress.com

subjID = fsss_subjID(JOB.subjID);
dir0=pwd;
NUMBER_OF_ROWS=5;
for n=1:numel(subjID)
subjid = subjID{n};
cd(fullfile(JOB.dir_base,subjid));
JOB.fname_png=[JOB.dir_png,'/',subjid,'_timecourse_',JOB.name_proc{:},'.png'];

NumProc = numel(JOB.fname_proc);
P=cell(1,NumProc);
y=cell(1,NumProc);
p=1;
P{p} = spm_vol(JOB.fname_proc{p});
xyz=[round(P{1}(1).dim(1)/2), round(P{1}(1).dim(2)/2), round(P{1}(1).dim(3)/(2))]';
y{p} = spm_get_data(P{p},xyz);
for p=2:NumProc
P{p} = spm_vol(JOB.fname_proc{p});
y{p} = spm_get_data(P{p},xyz);
end
dim=[P{1}(1).dim numel(P{1})];
numtp = min(100,dim(4));
if dim(4)>=200
tp=[round(dim(4)/4):round(dim(4)/4)+(numtp-1)];
else
tp=1:numtp;
end

ProcNames = JOB.name_proc;
figure;
set(gcf,'position',[841   206  1075  970]);
subplot(NUMBER_OF_ROWS,NumProc,[1:NumProc]);
yy=[];
for p=1:NumProc
yy=[yy zscore(y{p}(tp))];
end
plot(yy);
legend(ProcNames,'location','EastOutside')
title([subjid,':[',...
num2str(xyz(1)),',',num2str(xyz(2)),',',num2str(xyz(3)),'] vox'],...
'fontsize',16);
xlabel('TR'); ylabel('Z-score');

for p=1:NumProc
subplot(NUMBER_OF_ROWS,NumProc,NumProc+p)
x=spm_read_vols(P{p}(1));
imagesc(x(:,:,xyz(3))'); set(gca,'ydir','nor'); axis image; colormap(bone);
hold on; scatter(xyz(1),xyz(2),'r+');
axis off; title(ProcNames{p})
set(gca,'xtick',tp, 'xticklabel',tp);
end

% take 3 principle components
for p=1:NumProc
Y=zeros(dim(4),prod(dim(1:3)));
for j=tp
Y(j,:)=reshape(spm_read_vols(P{p}(j)),[1,prod(dim(1:3))]);
end
[~,score] = pca(Y);
PC{p} = score(:,[1:3]);
end

for k=1:3
subplot(NUMBER_OF_ROWS,NumProc,[1:NumProc]+(2+k-1)*NumProc);
ppc=[];
for p=1:NumProc
ppc=[ppc, PC{p}(:,k)];
end
plot(ppc);
xlabel('TR'); ylabel(['PC',num2str(k)]);
legend(ProcNames,'location','EastOutside')
end
xlim([tp(1) tp(end)])

screen2png(JOB.fname_png,120);
end
cd(dir0);
end
