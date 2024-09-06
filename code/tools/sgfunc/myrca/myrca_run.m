function job = myrca_run(job)
assert(exist('rcaRun.m','file'), ...
  'install (https://github.com/dmochow/rca) and set PATH for RCA toolbox')

% set output filenames:
nsubs = numel(job.subID);
thisjob = job;
thisjob.subid = job.subID{1};
thisjob = rmfield(thisjob,'rca');
thisjob.returnsuffixonly = true;
[~,~,thisjob] = feval(thisjob.prepfunc, thisjob);
fn_temp = thisjob.fn_data;

suffix_cond = cell2mat(job.conds);
suffix_rca = ['_rca',num2str(job.rca),'_',suffix_cond];

% check if all done:
alldone = true;
for isub = 1:nsubs
  fname_out{isub} = strrep(strrep(fn_temp, thisjob.subid, job.subID{isub}),...
    '.set',[suffix_rca,'.mat']);
  alldone = alldone && isfile(fname_out{isub});
end
if alldone
  return
end
% if not, let's do that:

% read epochs for EACH condition
data = {};
for isub = 1:nsubs
  
  % get RESP (false alarm rejected [optional], normalized) of EACH cond
  thisjob = job;
  thisjob.subid = job.subID{isub};
  thisjob = rmfield(thisjob,'rca');
%   thisjob = rmfield(thisjob,'conds');
  thisjob.returnresponly = true;
  [~, RESP, thisjob] = feval(thisjob.prepfunc, thisjob);
  assert(~sum(sum(cell2mat(cellfun(@isnan, RESP, 'uni',0)))))
  
  stims_trl = thisjob.eventlabels;
  stims_unq = unique(stims_trl);
  
  for icnd = 1:numel(stims_unq)
    idx_trl = find(ismember(stims_trl, stims_unq{icnd}));
    data{icnd,isub} = [];
    for itrl = idx_trl
      data{icnd,isub} = cat(3, data{icnd,isub}, RESP{itrl});
    end
  end
end


% run RCA
delete(gcp('nocreate'))
[data,W,A] = rcaRun(data, [], job.rca, [], [], 1, thisjob.chanlocs);
delete(gcp('nocreate'))
% save figure:
set(gcf,'colormap',parula,'visible','off')
dn_fig = fullfile(thisjob.dn_trfs,'figure','rca');
[~,~] = mkdir(dn_fig);
[~,suffix] = fileparts(thisjob.dn_trf);
suffix = [suffix, suffix_rca];
exportgraphics(gcf,fullfile(dn_fig,[suffix,'.png']))
close(gcf)


% save results as MAT files that you can use for MYTRF_RUN
for isub = 1:nsubs
  rca = data(:,isub);
  info = struct('stims_unq',{stims_unq'}, 'times',thisjob.times, ...
    'srate',thisjob.srate, 'chanlocs',thisjob.chanlocs, ...
    'W',W, 'A',A, 'job',job);
  save(fname_out{isub}, 'rca','info')
end


end
