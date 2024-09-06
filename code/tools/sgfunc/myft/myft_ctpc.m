function myft_ctpc(CFG)
% myft_ctpc(timefreq, CFG)
% (cc) sgKIM, 2019.

CFG.savemat = 0;
CFG.precision = 'single';
[timefreq, CFG] = myft_freqanalysis(CFG);

iCond = 1;
if ischar(timefreq.(CFG.cond(iCond).name))
  load(timefreq.(CFG.cond(iCond).name),'tfr')
else
  tfr = timefreq.(CFG.cond(iCond).name);
end
d = size(tfr.fourierspctrm);
% #trials x #chans x #freqs x #times

% CTPC-within:
tic
ctpc_w = zeros(numel(CFG.cond), d(2), d(3), d(4));
nTrl = zeros(1, numel(CFG.cond));
for iCond = 1:numel(CFG.cond)
  if ischar(timefreq.(CFG.cond(iCond).name))
    load(timefreq.(CFG.cond(iCond).name),'tfr')
  else
    tfr = timefreq.(CFG.cond(iCond).name);
  end
  F = tfr.fourierspctrm;
  nTrl(iCond) = size(F,1);
  ctpc_w(iCond,:,:,:) = myft_compute_ctpc(F);
end
toc

%{
NOTE: without matching the # of trials, that is, if I use all the trials
together than the CTPC is estimated much weaker than when the # of trials
was matched. Because 30-50 trials are still small number for stable
estimation of phase-locking, so in this case it would be better to
"control" the # of trials, too. (although it would take x50 longer)
%}
% CTPC-across-exemplars
if isfield(CFG,'ctpc_x_conds')
  tic
  ctpc_x = zeros(max(CFG.ctpc_x_conds), d(2), d(3), d(4));
  % across exemplars belong to the same condition
  for jXCond = 1:max(CFG.ctpc_x_conds)
    idxCond = find(CFG.ctpc_x_conds == jXCond);
    nTrl_per_cond = round( mean(nTrl(idxCond))/numel(idxCond) );
    % randomly sample the equal # of trials from conditions:
    ctpc_x_rand = zeros(100, d(2), d(3), d(4));
    for iRand = 1:100
      F = [];
      for j = 1:numel(idxCond)
        jCond = idxCond(j);
        if ischar(timefreq.(CFG.cond(jCond).name))
          load(timefreq.(CFG.cond(jCond).name),'tfr')
        else
          tfr = timefreq.(CFG.cond(jCond).name);
        end
        f = tfr.fourierspctrm;
        idxTrl = randperm(size(f,1));
        F = cat(1, F, f(idxTrl(1:nTrl_per_cond),:,:,:));
      end
      ctpc_x_rand(iRand,:,:,:) = myft_compute_ctpc(F);
    end
    ctpc_x(jXCond,:,:,:) = squeeze(mean(ctpc_x_rand,1)); % average
  end
  toc
end

tfr = rmfield(tfr, 'fourierspctrm'); % only meta info
fn_mat = [CFG.prefix,'_ctpc.mat'];
save(fn_mat,'tfr','CFG','ctpc_w')
if exist('ctpc_x','var')
  save(fn_mat,'-append','ctpc_x')
end
ls(fn_mat)
end