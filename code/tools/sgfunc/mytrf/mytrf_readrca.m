function [RESP, job] = mytrf_readrca(job)

thisjob = job;
thisjob = rmfield(thisjob,'rca');
thisjob.returnsuffixonly = true;
[~,thisjob] = mytrf_prepresp(thisjob);

suffix_cond = cell2mat(job.conds);
suffix_rca = ['_rca',num2str(job.rca),'_',suffix_cond];

job.fn_data = strrep(thisjob.fn_data, '.set',[suffix_rca,'.mat']);
load(job.fn_data, 'rca','info')

% find labels for given condition codes:
if isfield(job,'conds')
  include = false(size(info.stims_unq));
  for icnd = 1:numel(job.conds)
    include = include | contains(info.stims_unq, job.conds{icnd});
  end
  if ~any(include)
    error('No trials found for given job.conds!')
  end
  % select for given conditions:
  info.stims_unq = info.stims_unq(include);
  rca = rca(include,1);
end


% reshape into 1-D cell array:
ntrials = cell2mat(cellfun(@(x) size(x,3), rca, 'uni',0));
job.trials = sum(ntrials);
RESP = cell(job.trials,1);
job.eventlabels = cell(job.trials,1);
k = 0;
for icnd = 1:numel(rca)
  for itrl = 1:size(rca{icnd},3)
    RESP{itrl+k,1} = rca{icnd}(:,:,itrl);
    job.eventlabels{itrl+k,1} = info.stims_unq{icnd};
  end
  k = k + ntrials(icnd);
end
job.times = info.times;
job.srate = info.srate;
job.dn_trf = fileparts(job.fn_data);
job.chanlocs = info.chanlocs;
job.W = info.W;
job.A = info.A;

end