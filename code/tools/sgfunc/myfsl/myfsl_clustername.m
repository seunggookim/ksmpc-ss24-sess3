function JOB = myfsl_clustername(JOB)
% JOB = myfsl_clustername(JOB)
%
% JOB requires:
%  .fname_pval [!!!] must be in MNI152 space [!!!]
% (.alpha)
%
% (cc) 2015, sgKIM.   solleo@gmail.com

[dir1,name1,ext1] = fileparts_gz(JOB.fname_pval);
if ~isfield(JOB,'alpha'), JOB.alpha=0.05; end
thres = num2str(1-JOB.alpha);
JOB.fname_txt    = fullfile(dir1,[name1,'_cluster_thrs_',thres,'.txt']);
JOB.fname_clusid = fullfile(dir1,[name1,'_clusterid_thrs_',thres,'.nii']);
myunix(['FSLOUTPUTTYPE=NIFTI_GZ; ', ...
' cluster -i ',JOB.fname_pval,' --thresh=',thres,' --mm >',JOB.fname_txt,'; ', ...
' cluster -i ',JOB.fname_pval,' --thresh=',thres,' --mm --oindex=',JOB.fname_clusid]);
M = dlmread(JOB.fname_txt,'\t',1,0);
for j=1:max(M(:,1))
JOB.clus(j).id = M(j,1);
JOB.clus(j).numVox = M(j,2);
JOB.clus(j).pval = 1-M(j,3);

strct = myfsl_atlasquery(M(j,4:6));
JOB.clus(j).peak_coord_mm = M(j,4:6);
JOB.clus(j).peak_name = strct.name;
JOB.clus(j).peak_prob_prcnt = strct.prob;

strct = myfsl_atlasquery(M(j,7:9));
JOB.clus(j).cog_coord_mm  = M(j,7:9);
JOB.clus(j).cog_name  = strct.name;
JOB.clus(j).cog_prob_prcnt = strct.prob;
end
JOB.fname_mat = fullfile(dir1,[name1,'_cluster_thrs_',thres,'.mat']);
clus = JOB.clus;
save(JOB.fname_mat,'clus');

end
