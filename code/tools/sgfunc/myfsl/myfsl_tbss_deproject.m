function myfsl_tbss_deproject (JOB)
% myfsl_tbss_deproject (JOB)
%
% (cc) 2015, sgKIM.

if ~isfield(JOB,'alpha'), JOB.alpha=0.05; end
thres = num2str(1-JOB.alpha);

% 1. take p-value image and create cluster-id maps
JOB.fname_pval = fullfile(JOB.dir_glm,JOB.name_pval);
JOB = myfsl_clustername(JOB);

% 2. run tbss_deproject (takes time to create inverse warp files (why???))
myunix(['cd ',JOB.dir_stats,'; tbss_deproject ',JOB.fname_clusid,' 2 -n']);


end
