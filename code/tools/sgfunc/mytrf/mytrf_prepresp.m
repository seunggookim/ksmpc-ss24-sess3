function [RESP, job] = mytrf_prepresp(job)
% [RESP, job] = mytrf_prepresp(job)
%

%% Option
if isfield(job,'reref')
  suffix_reref = ['_',job.reref];
else
  suffix_reref = '';
end
if ~isfield(job,'rejresp')
  job.rejresp = 'none';
end
if ~isfield(job,'normresp')
  job.normresp = true;
end
if ~isfield(job,'autoreject')
  job.autoreject = false;
end

%% Load data & finish the preproc:
job.fn_data = fullfile(job.dn_trf, sprintf(...
  '%s_ds%gHz_bpf%g-%gHz_toi%g-%g%s_normresp%g_autoreject%g/epochs.set', ...
  job.suffix_data, job.ds_Hz, job.bpf_Hz, job.toi_s, suffix_reref, ...
  job.normresp, job.autoreject));
if isfield(job,'returnsuffixonly') && job.returnsuffixonly
  RESP=[];
  return
end
[p1,f1,e1] = fileparts(job.fn_data);
[~,~] = mkdir(p1);
fn_log = fullfile(p1,[f1,'.log']);

% because triggers at t=0 are needed for merging and re-epoching:
toi = [min(job.toi_s(1),-1) job.toi_s(2)];

if isfile(job.fn_data)
  eeg = pop_loadset(job.fn_data);
else
  if isfile(fn_log); delete(fn_log); end
  diary(fn_log)
  fprintf('[%s]%s\n',datestr(now,31), repmat('-',[1 50]))
  disp(job)
  fprintf('[%s] DATA MERGING/PREPROCESSING starts.\n',datestr(now,31))
  t0 = tic;
  
  % Merge all runs for all conditions:
  for ifile = 1:numel(job.fnames)
    % read one run:
    fn = fullfile(job.fnames{ifile});
    eeg1 = pop_loadset(fn);
    
    % down-sampling:
    eeg1 = pop_resample(eeg1, job.ds_Hz);
    
    % filtering:
    if isfield(job,'bpf_Hz')
      if job.bpf_Hz(1)==0 && ~isinf(job.bpf_Hz(2))
        eeg1 = pop_eegfiltnew(eeg1, [], job.bpf_Hz(2));
      elseif job.bpf_Hz(1)>0 && isinf(job.bpf_Hz(2))
        eeg1 = pop_eegfiltnew(eeg1, job.bpf_Hz(1), []);
      else
        eeg1 = pop_eegfiltnew(eeg1, job.bpf_Hz(1), job.bpf_Hz(2));
      end
    end
    
    % epoching:
    idx_events = find(ismember({eeg1.event.code},'Stimulus'));
    fprintf('%i "Stimulus" events found\n',numel(idx_events))
    eeg1 = pop_epoch(eeg1, [], toi, 'eventindices', idx_events);
    
    % merge:
    if ~exist('eeg','var')
      eeg = eeg1;
    else
      eeg = pop_mergeset(eeg, eeg1, 1);
    end
  end
  clear eeg1
  
  % auto-rejection of bad trials and interpolation of bad channels
  if isfield(job,'autoreject') && job.autoreject
    suffix = sprintf('%s_ds%gHz_bpf%g-%gHz', ...
      job.suffix_data, job.ds_Hz, job.bpf_Hz);
    eeg = myeeg_autoreject(eeg, true, suffix);
  end
  
  % re-referencing:
  if isfield(job,'reref')
    switch (job.reref)
      case 'avgmastoid' % this works only when ref=RM
        thischan = ismember({eeg.chanlocs.labels},'LM');
        eeg.data = eeg.data - 0.5*eeg.data(thischan,:,:);
        eeg.ref = 'avgmastoid';
      case 'common'
        eeg = pop_reref(eeg, []);
      case 'laplacian'
        eeg1 = eeg;
        for itrl = 1:eeg.trials
          eeg1.data = eeg.data(:,:,itrl);
          eeg.data(:,:,itrl) = eeg_laplac(eeg1, 1);
        end
        eeg.ref = 'laplacian';
      otherwise
        error('UNKNONW OPTION FOR JOB.REREF = %s', job.reref)
    end
  end
  
  eeg.filepath = p1;
  eeg.filename = [f1,e1];
  eeg.setname = 'ready for mTRF';
  pop_saveset(eeg, 'savemode','resave');
  fprintf('[%s] DATA MERGING/PREPROCESSING ends. ',datestr(now,31))
  toc(t0)
  diary off
end

if isfield(job,'saveepochsonly') && job.saveepochsonly
  RESP = [];
  return
end

% AFTER REJECTION, we need to rely on EEG.EPOCH.EVENT*
if iscell(eeg.epoch(1).eventtype)
  epochevents = [eeg.epoch.eventtype];
  epochcodes = [eeg.epoch.eventcode];
else
  epochevents = {eeg.epoch.eventtype};
  epochcodes = {eeg.epoch.eventcode};
end
epochstims = epochevents(ismember(epochcodes,'Stimulus'));
assert(numel(epochstims) == eeg.trials, '#epochstims =/= #trials?')

% find labels for given condition codes:
if isfield(job,'conds')
  include = false(size(epochstims));
  for icond = 1:numel(job.conds)
    include = include | contains(epochstims, job.conds{icond});
  end
  if ~any(include)
    error('No trials found for given job.conds!')
  end
  % select for given conditions:
  eeg = pop_select(eeg, 'trial', find(include));
end

% check for "Response" based on EEG.EPOCH.EVENT*
if iscell(eeg.epoch(1).eventcode)
  idx_trial = find(cell2mat(cellfun(@(x) any(contains(x, 'Response')), ...
    {eeg.epoch.eventcode},'uni',0)));
else
  idx_trial = find(ismember({eeg.epoch.eventcode}, 'Response'));
end

if ~isempty(idx_trial)
  diary(fn_log)
  fprintf('[%s]%s\n',datestr(now,31), repmat('-',[1 50]))
  disp(job)
  warning('"Response" found in %i/%i trials',numel(idx_trial), eeg.trials);
  inout = {'OUT','IN'};
  isin = false(1,numel(idx_trial));
  for i = 1:numel(idx_trial)
    idx = idx_trial(i);
    idx_resp = find(ismember(eeg.epoch(idx).eventcode,'Response'));
    for j = idx_resp
      t = eeg.epoch(idx).eventlatency{j};
      isin(i) = isin(i) || (t/1000>=job.toi_s(1) && t/1000<=job.toi_s(2));
      fprintf('[epoch%i:%gms] "%s" (%s) \n', ...
        idx, t, eeg.epoch(idx).eventtype{j}, inout{isin(i)+1})
    end
  end
  switch lower(job.rejresp)
    case {'onlyin','inonly','onlywithin','onlywithintoi'}
      if sum(isin)
        fprintf(['REJECTING ONLY trials with responses WITHIN',...
          '[%g,%g] sec (%i/%i event[s]): trial#'],...
          job.toi_s, sum(isin), numel(isin));
        trials2rej = idx_trial(isin);
      else
        fprintf('REJECTING NOTHING\n')
        trials2rej = [];
      end
      
    case {'all','everything','inandout','in&out'}
      fprintf(['REJECTING ANY trials with responses WITHIN ',...
        '[-1,%g] sec (%i event[s]): trial#'],...
        job.toi_s(2), numel(isin));
      trials2rej = idx_trial;
      
    otherwise
      fprintf('REJECTING NOTHING for false alarm responses\n')
      trials2rej = [];
  end
  if ~isempty(trials2rej)
    fprintf('%i, ',trials2rej);
    fprintf('\b\b\n');
    [eeg] = pop_select(eeg, 'notrial',trials2rej);
  end
  diary off
end

% AFTER REJECTION, we need to rely on EEG.EPOCH.EVENT*
if iscell(eeg.epoch(1).eventtype)
  epochevents = [eeg.epoch.eventtype];
  epochcodes = [eeg.epoch.eventcode];
else
  epochevents = {eeg.epoch.eventtype};
  epochcodes = {eeg.epoch.eventcode};
end
epochstims = epochevents(ismember(epochcodes,'Stimulus'));
assert(numel(epochstims) == eeg.trials, '#epochstims =/= #trials?')

% re-epoching for toi_s(1) > 0
[eeg] = pop_select(eeg, 'time', job.toi_s);

% Get metainfo fields:
job.trials = eeg.trials;
job.chanlocs = eeg.chanlocs;
job.times = eeg.times;
job.eventlabels = epochstims;
job.srate = eeg.srate;
% job.fn_data = fn_data;
job.dn_trf = fileparts(job.fn_data);

% Return RESP
if job.normresp % Normalization of EEG responses with global mean/std:
  globalmean = mean(reshape(eeg.data, ...
    [eeg.nbchan, eeg.pnts*eeg.trials]),2);
  globalstd = std(reshape(eeg.data, ...
    [eeg.nbchan, eeg.pnts*eeg.trials]),[],2);
  RESP = cellfun(@(x) ...
    ((squeeze(eeg.data(:,:,x)-globalmean))./globalstd)', ...
    num2cell(1:eeg.trials)', 'uni',0);
else
  RESP = cellfun(@(x) squeeze(eeg.data(:,:,x))', ...
    num2cell(1:eeg.trials)', 'uni',0);
end

end
