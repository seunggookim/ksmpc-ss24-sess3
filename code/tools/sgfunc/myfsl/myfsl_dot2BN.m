function myfsl_dot2BN
%% Read conn matrix
if exist('fdt_matrix3.dot','file') && ~exist('fdt_matrix3.mat','file')
 dot=load('fdt_matrix3.dot','-ascii');
 X=spconvert(dot);
 save fdt_matrix3.mat X
 clear dot X
 delete('fdt_matrix3.dot');
end
load('fdt_matrix3.mat')
degree = full(sum(X));

%%
% nii=load_untouch_nii('all.nii.gz');
coords=load('coords_for_fdt_matrix3','-ascii');
% D=nii.img*0;
% for a=1:size(coords,1)
%  i=coords(a,1)+1;
%  j=coords(a,2)+1;
%  k=coords(a,3)+1;
%  D(i,j,k)=degree(a);
% end
% not very informative...
%%
if ~exist('../data/daparc+aseg.BN_Atlas.nii','file')
 if ~exist('../../fs6/mri/aparc+aseg.BN_Atlas.nii','file')
  unix(['mri_convert ../../fs6/mri/aparc+aseg.BN_Atlas.mgz ', ...
   '../../fs6/mri/aparc+aseg.BN_Atlas.nii']);
 end
 if ~exist('../../fs6/mri/brain.nii','file')
  unix(['mri_convert ../../fs6/mri/bain.mgz ../../fs6/mri/brain.nii']);
 end
 job1=[];
 job1.prefix='d';
 job1.interp=0;
 job1.fname_moving='../../fs6/mri/brain.nii';
 job1.fname_fixed ='../../fdt/data/bmnodif.nii';
 job1.fname_others='../../fs6/mri/aparc+aseg.BN_Atlas.nii';
 myspm_coreg(job1)
 unix(['mv ../../fs6/mri/daparc+aseg.BN_Atlas.nii ', ...
  '../data/daparc+aseg.BN_Atlas.nii']);
end
%%
BN=load_untouch_nii('../data/daparc+aseg.BN_Atlas.nii');
idx_BN=zeros(size(coords,1),1);
for a=1:size(coords,1)
 i=coords(a,1)+1;
 j=coords(a,2)+1;
 k=coords(a,3)+1;
 idx_BN(a)=BN.img(i,j,k);
end
disp(...
 ['BN voxels without tracks= ', ...
 num2str(sum(~degree)),' / ',num2str(size(coords,1)), ...
 ' (',num2str(sum(~degree)/size(coords,1)*100),' %)' ]);
%%
T=readtable('~/Dropbox/BN_Atlas/BN_Atlas_246_COT.xlsx');
M=zeros(246);
for i=1:245
 for j=(i+1):246
  idx_i = idx_BN==T.Label(i);
  idx_j = idx_BN==T.Label(j);
  M(i,j) = full(sum(sum(X(idx_i,idx_j)))) / (sum(idx_i)*sum(idx_j));
  M(j,i) = M(i,j);
 end
end
save('../../BN_vol/FDT.mat','M')

% cfg=[];
% %cfg.colorbartitle=['PCC_WCC',num2str(j),'_',subj];
% %cfg.fname_png=['PCC_cor/scale_',num2str(j),'_',subj,'.png'];
% cfg.caxis=[0 500];
% fsss_view_BN(nansum(M), cfg);

%%
end