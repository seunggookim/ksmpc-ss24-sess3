function [CV,job] = mytrf_cv(STIM, RESP, job)
%mytrf_cv - a wrapper for dual-loop CV using mTRFcrossval()
% [CV,job] = mytrf_cv(STIM, RESP, job)
%
%[intput]
% STIM          {nx1} a cell array of stimuli for n trials
% RESP          {nx1} a cell array of responses for n trials
% job           (1x1) a job configuration structure
%  .nfolds      [1x1] (default = 10)
%  .outercv     '1xc' 'stim' (default) | 'trial'
% (.band)       [1xb] # of regressors per band; eg. [1,2] for [a1,b1,b2]
%  .optmethod   '1xc' 'minmse' (default) | 'maxr'
%  .zeropad     [1x1] 0 (default) | 1
%  .lag_ms      
%  .srate
%  .lambdas
% (.debug)
% (.dn_trf)
%
%[output]
% CV             (1x1) an output structure
%  .mse_opt      MSE from optimization (inner-CV)
%  .mse_udf      "underfitting MSE" (when prediction is a flat zero line)
%  .r_opt        r from optimization (inner-CV)
%  .lambda_opt   optimized lambda
% (.lambdaset)   optimized lambdaset (for banded ridge)
%  .model        models from evaluation (outer-CV)
%  .mse_eva      MSE from evaluation (outer-CV)
%  .r_eva        r from evaluation (outer-CV)
%
% (cc) 2021, sgKIM.


%% Set up parameters:
if ~exist('job','var'), job = []; end
if ~isfield(job,'nfolds'), job.nfolds = 10; end
if ~isfield(job,'outercv'), job.outercv = 'stim'; end
if ~isfield(job,'optmethod'), job.optmethod = 'minmse'; end
if ~isfield(job,'zeropad'), job.zeropad = 0; end
if ~isfield(job,'lag_ms'), job.lag_ms = [-100 400]; end
if ~isfield(job,'lambdas'), job.lambdas = 10.^(-10:5); end
if ~isfield(job,'dir'), job.dir = 1; end
if ~isfield(job,'fast')
  if job.dir == 1, job.fast = 1; else, job.fast = 0; end
end
if ~isfield(job,'type')
  if job.dir == 1, job.type = 'multi'; else, job.type = 'single'; end
end
if ~isfield(job,'fitmethod'), job.fitmethod = 'ridge'; end
if job.dir == -1, dirtxt = 'dec'; elseif job.dir == 1, dirtxt = 'enc'; end

%% DEFINE output suffix
suffix = '';
if isfield(job,'suffix_trf'), suffix = [suffix,'_',job.suffix_trf]; end
if isfield(job,'randstim') && ~strcmp(job.randstim,'none')
  suffix = [suffix,'_rnd-',job.randstim];
end
if isfield(job,'conds') && ~isfield(job,'stimcond')
  job.stimcond = cell2mat(job.conds);
end
if isfield(job,'stimcond')
  suffix = [suffix,'_cnd-',job.stimcond];
else
  suffix = [suffix,'_cnd-all'];
end
if isfield(job,'band')
  suffix = [suffix,'_bnd-',num2str(numel(job.band))];
end
if ~strcmp(job.fitmethod,'ridge')
  suffix = [suffix,'_fit-',job.fitmethod];
end
if isfield(job,'rejresp')
  suffix = [suffix,'_rejresp-',job.rejresp];
else
  suffix = [suffix,'_rejresp-none'];
end
  
if isfield(job,'dn_trf')
  job.dn_out = fullfile(job.dn_trf, sprintf(...
    '%s_ocv-%s_opt-%s_zp%g_lag%g-%gms_ltn%g-%g%s',...
    dirtxt, job.outercv, job.optmethod, job.zeropad, job.lag_ms, ...
    log10(job.lambdas([1 end])), suffix ));
  [~,~] = mkdir(job.dn_out);
  job.fn_cv = fullfile(job.dn_out,'cv.mat');
  job.fn_log = fullfile(job.dn_out,'cv.log');
  if isfile(job.fn_cv)
    load(job.fn_cv,'CV','job')
    return % if it's already done
  end
end
%%
diary(job.fn_log)
fprintf('[%s]%s\n',datestr(now,31), repmat('-',[1 50]))
disp(job)

%%
mse_tr = []; mse_opt = []; mse_eva = []; 
r_tr   = []; r_opt   = []; r_eva   = [];
lambda_opt = []; mdls = [];

% Set up outer-CV partitions[*]:
if ~isfield(job,'eventlabels_outercv')
  job.eventlabels_outercv = job.eventlabels;
end
switch job.outercv
  case 'stim'
    % k-fold random partition of stimuli
    % (eg. stim1 could be only either in a training set or a test set)
    uniqueevents = unique(job.eventlabels_outercv);
    cvpart = cvpartition(numel(uniqueevents),'kfold',job.nfolds);
  case 'trial'
    % k-fold random partition of trials
    % (eg. stim1 could be both in a training set and a test set)
    cvpart = cvpartition(numel(STIM),'kfold',job.nfolds);
  otherwise
    error('%s UNKNOWN',job.outercv)
end

for icv = 1:job.nfolds % OUTER-CV (k-fold)
  % partition inner-CV and holdout (testing) data
  switch job.outercv
    case 'stim'
      events_tr = uniqueevents(training(cvpart,icv));
      ind_tr = ismember(job.eventlabels_outercv, events_tr);
      events_te = uniqueevents(test(cvpart,icv));
      ind_te = ismember(job.eventlabels_outercv, events_te);
    case 'trial'
      ind_tr = training(cvpart, icv);
      ind_te = test(cvpart, icv);
  end
  
  % optimize on validation data (INNER-CV: LOOCV)
  fprintf('[OUTER-CV:%i/%i] Optimizing...', icv, job.nfolds)
  t1=tic;
  if isfield(job,'band')
    stat_opt = mTRFcrossval(STIM(ind_tr),RESP(ind_tr),job.srate,job.dir,...
      job.lag_ms(1),job.lag_ms(2),job.lambdas, 'zeropad',job.zeropad,...
      'band',job.band, 'method',job.fitmethod, 'verbose',0, ...
      'fast',job.fast, 'type',job.type);
  else
    stat_opt = mTRFcrossval(STIM(ind_tr),RESP(ind_tr),job.srate,job.dir,...
      job.lag_ms(1),job.lag_ms(2),job.lambdas,'zeropad',job.zeropad,...
      'method',job.fitmethod, 'verbose',0, 'fast',job.fast, ...
      'type',job.type);
  end
  fprintf('took %i sec.\n', round(toc(t1)));
  
  mse_opt(icv,:) = mean(mean(stat_opt.err,1),3);  % average across channels
  mse_udf(icv,:) = mean(cell2mat(cellfun(...      % upper bound of MSE
    @(x) mean(mean((x-mean(x)).^2,1),2), RESP(ind_tr),'UniformOutput',0)));
  r_opt(icv,:) = mean(mean(stat_opt.r,1),3);      % average across channels
  switch job.optmethod
    case 'minmse'
      [~,idx] = min(mse_opt(icv,:));
    case 'maxr'
      [~,idx] = max(r_opt(icv,:));
  end
  if isfield(job,'band')
    lambda_opt(icv,:) = stat_opt.lambdaset(idx,:);
  else
    lambda_opt(icv,1) = job.lambdas(idx);
  end
  
  % fit on training/validation data [+]
  if isfield(job,'band')
    mdl = mTRFtrain(STIM(ind_tr), RESP(ind_tr), job.srate,1,...
      job.lag_ms(1),job.lag_ms(2), lambda_opt(icv,:), ...
      'zeropad',job.zeropad, 'band',job.band, 'method',job.fitmethod, ...
      'verbose',0);
  else
    mdl = mTRFtrain(STIM(ind_tr), RESP(ind_tr), job.srate,1,...
      job.lag_ms(1),job.lag_ms(2), lambda_opt(icv), ...
      'zeropad',job.zeropad, 'method',job.fitmethod, 'verbose',0);
  end
  mdls = [mdls mdl];
  
  % get training error
  [~,stat_tr] = mTRFpredict(STIM(ind_tr), RESP(ind_tr), mdl, ...
    'zeropad',job.zeropad, 'verbose',0);
  mse_tr(icv,:) = mean(stat_tr.err,1);     % average across trials
  r_tr(icv,:) = mean(stat_tr.r,1);         % average across trials
  fprintf('training r = %.4f, error = %.4f\n', ...
    mean(r_tr(icv,:)), mean(mse_tr(icv,:)));
  
  % evaluate on testing (holdout) data
  [~,stat_eval] = mTRFpredict(STIM(ind_te), RESP(ind_te), mdl, ...
    'zeropad',job.zeropad, 'verbose',0);
  mse_eva(icv,:) = mean(stat_eval.err,1);  % average across trials
  r_eva(icv,:) = mean(stat_eval.r,1);      % average across trials
  fprintf('testing r = %.4f, error = %.4f\n', ...
    mean(r_eva(icv,:)), mean(mse_eva(icv,:)));
  
end


%% return everything:
CV = [];
CV.mse_opt = mse_opt;
CV.lambda_opt = lambda_opt;
if isfield(stat_opt,'lambdaset')
  CV.lambdaset = stat_opt.lambdaset;
end
CV.mse_udf = mse_udf;
CV.r_opt = r_opt;
CV.model = mdls;
CV.mse_tr = mse_tr;
CV.r_tr = r_tr;
CV.mse_eva = mse_eva;
CV.r_eva = r_eva;

if isfield(job,'fn_cv')
  save(job.fn_cv, 'CV','job')
end
diary off

if ~(isfield(job,'debug') && job.debug)
  return
end


%% [DEBUGGING MODE] detect underfitting & overfitting:
%{
usually: MSE_tr < MSE_opt < MSE_eva < MSE_udf
(1) Everything must be larger than MSE_udf: so MSE_eva <= MSE_udf
(2) MSE_eva < MSE_tr: possible with regularization
%}

figure
for icv = 1:job.nfolds
  mse_opt = CV.mse_opt(icv,:);
  mse_eva = mean(CV.mse_eva(icv,:),2);
  mse_tr = mean(CV.mse_tr(icv,:),2);
  mse_udf = CV.mse_udf(icv);
  
  p = nextpow2(job.nfolds);
  subplot(p,p,icv)
  plot(job.lambdas, mse_opt,'k')
  hold on
  ind = job.lambdas == CV.lambda_opt(icv);
  h(1) = scatter(CV.lambda_opt(icv), mse_opt(ind),'ro');
  h(2) = scatter(CV.lambda_opt(icv), mse_tr,'gx');
  h(3) = scatter(CV.lambda_opt(icv), mse_eva,'bo');
  h(4) = scatter(CV.lambda_opt(icv), mse_udf,'cx');
  set(gca,'xscale','log')
  axis square
  
  title(['cv',num2str(icv)])
  assert(mse_eva<=mse_udf, 'mse_eval > mse_udf?')  
end
legend(h, {'opt','tr','eval','udf'})


end


%{
MODEL AVERAGING OR THIS? I THINK MODEL AVERAGING IS OKAY FOR LINEAR MODELS
% if isfield(job,'band')
%   % train the whole data with mean lambda:
%   model = mTRFtrain(STIM,RESP,job.srate,1,...
%     job.lag_ms(1),job.lag_ms(2), geomean(lambda_opt,1),...
%     'zeropad',job.zeropad,'band',job.band);
% else
%   % train the whole data with mean lambda:
%   model = mTRFtrain(STIM,RESP,job.srate,1,...
%     job.lag_ms(1),job.lag_ms(2), geomean(lambda_opt,1),...
%     'zeropad',job.zeropad);
% end
%(or just model averaging?)

%}

%{ 
[*] OUTER-CV
Stratified CV: this is to minimize the influence from the bias in samples.
For a binary classifier, if training sets have more of class A and less of
class B, then the classifier could be biased towards class A (really?),
thus performing non-optimally. To avoid this, one could make sure that each
of folds have equal proportions of all classes (ie. pseudorandomization).
For regression, a similarity measure of true values is used?
%}

%{
[+] re-fitting to the training/validation set
: This could be even more overfitting to the validation set IF we
partitioned the data only once. But here we use LOOCV, so each trial was
used as a validation set once. And the optimal lambda was deteremined from
the across-inner-fold-averaged error curve (which is wrong?).


VALIDATION SET (Bishop, 1995)
"""
The performance of the networks is then compared by evaluating the error
function using an independent validation set, and the network having the
smallest error with respect to the validation set is selected. This
approach is called the hold out method. Since this procedure can itself
lead to some overfitting to the validation set, the performance of the
selected network should be confirmed by measuring its performance on a
third independent set of data called a test set.
"""
%}
