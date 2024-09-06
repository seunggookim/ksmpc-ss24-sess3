function myspm_mni_to_BN_vol(fname_epi,figs)
% myspm_mni_to_BN_vol(fname_epi)
% (cc) sgKIM.

if ~exist('figs','var'), figs=0; end
[dir1,f1,~]=myfileparts(fname_epi);
if ~exist('/tmp/BN1.nii.gz','file')
 unix('cp ~/Dropbox/BN_Atlas/BN_Atlas_246_1mm.nii.gz /tmp/BN1.nii.gz');
end
if ~exist('/tmp/BN1.nii','file')
 unix('gunzip /tmp/BN1.nii.gz');
end
if ~exist('/tmp/rBN1.nii','file')
 job1=[];
 job1.fname_ref=fname_epi;
 job1.fname_src='/tmp/BN1.nii';
 job1.interp=0;
 myspm_reslice(job1);
end
%T=readtable('~/Dropbox/BN_Atlas/BN_Atlas_246_COT.xlsx');

%% now averaging within each parcellation
nii=load_untouch_nii(fname_epi);
nii.img=double(nii.img);
maxT=size(nii.img,4);
Y=zeros(maxT,246);
BN=load_untouch_nii('/tmp/rBN1.nii');
y=img2y(nii.img);
for t=1:maxT
 for n=1:246
  ind=BN.img(:)==n;
  Y(t,n)=nanmean(y(t,ind));
 end
end
tSNR=mean(Y)./std(Y);
save([dir1,'/BN_',f1,'.mat'],'Y','tSNR')

%% Show how it's done:
if figs
 cfg=[];
 cfg.colorbartitle='tSNR';
 cfg.fname_png=[dir1,'/fig_',f1,'_tSNR.png'];
 if ~exist(cfg.fname_png,'file')
  fsss_view_BN(tSNR, cfg);
  I=imread(cfg.fname_png);
  I=imresize(I, 0.2);
  imwrite(I,cfg.fname_png);
 end
else
end
%%


end