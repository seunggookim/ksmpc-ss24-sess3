% function myspm_motion_orthogonality(JOB)

% dir2=JOB.dir_glm;

load([dir2,'SPM.mat']);
mus=SPM.xX.X(:,1:4);
mus=[mus mean(mus,2)];
gm=load([dir1,'cc_gm.txt'],'-ascii');
cc=load([dir1,'cc_n3d1v1b0.01-Inf_eigenvec.txt']);
rp=load([dir1,'rp_afunc8.txt']);

[R,P]=corr(mus,[gm cc rp]);
figure;
imagesc(R);axis image; colorbar
caxis([-.5 .5]); colormap bipolar
set(gca,'ytick',1:5,'yticklabel',{'FC','FD','BC','BD','mean'}, ...
 'xtick',1:4,'xticklabel',{'gm','cc1','cc2','cc3','rp..'});
% end