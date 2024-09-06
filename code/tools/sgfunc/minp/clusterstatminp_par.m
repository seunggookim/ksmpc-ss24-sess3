function [stat, cfg] = clusterstatminp_par(job, statobs, statrnd)
%clusterstatminp_par Performs cluster-based permutation tests with min(p)
%
%NOTES
% - CURRENTLY tested for: [1] 3-D volumetric or 2-D surface data,
% [2] one-sided or two-sided one-sample tests
% - THIS USES A CUSTOM VERSION of FieldTrip's CLUSTERSTAT(). Runs without a
% FieldTrip installation pathed.
%
%SYNTAX
% [stat, cfg] = clusterstatminp_par(job, statobs, statrnd)
%
%INPUTS
% job
%  .dim           [1*#dim|num]  # of voxels along x, y, z-directions
%  .inside        [#vox*1:logical]  true=in, false=out
%  .tail          [1*1|num]  1=positive, -1=negative, 0=both
%  .clusteralphas [1*#alphas|num]  clutser-alphas to run
%  .clusterconns  [1*#conn|num]  connectivity criterion
%        2D: 4 (edge), 8 (corner)
%        3D: 6 (surface), 18 (edge) or 26 (corner)
%  .connectivity  [#vox*#vox:sparse-logical]
%
% statobs [#vox*1|num] if this is a GpuArray it will run on GPU.
%
% statrnd [#vox*#rands|num]
%
%OUTPUTS
%
%
%REF: Geerligs, L., & Maris, E. (2021). Improving the sensitivity of
% cluster-based statistics for functional magnetic resonance imaging data.
% Human brain mapping, 42(9), 2746-2765. https://doi.org/10.1002/hbm.25399
%
% Implementing the "min(p)" method combining multiple cluster-defining
% thresholds (CDTs; controlling height) and types of connectivity
% (controlling cluster size).
%
% (cc) 2021, dr.seunggoo.kim@gmail.com

if isa(statobs,'gpuArray')
  arrayclass = 'gpuArray';
else
  arrayclass = 'double';
end

if ~isfield(job,'clusterstatistic')
  job.clusterstatistic = 'maxsum';
end

[nvox, nrnd] = size(statrnd);

posdistributionminp = ones(1,nrnd,arrayclass);
negdistributionminp = ones(1,nrnd,arrayclass);

posobsminp = nan(nvox,1,arrayclass);
negobsminp = nan(nvox,1,arrayclass);

% by default use MAXIMAL # of workers on the system
c = parcluster('local');
delete(gcp('nocreate'));
g = parpool(c, c.NumWorkers);  

for icdt = 1:numel(job.clusteralphas)
  for icon = 1:numel(job.clusterconns)
    cfg = struct( 'tail',job.tail, 'dim',job.dim, 'inside',job.inside, ...
      'connectivity',job.connectivity, ...
      'orderedstats','no', 'multivariate','no', 'minnbchan',0, ...
      'wcm_weight',1, 'clusterstatistic',job.clusterstatistic, ...
      'clusterthreshold','nonparametric_common',...
      'clusteralpha', job.clusteralphas(icdt), 'clustercritval',[], ...
      'clustertail',job.tail, 'clusterconn', job.clusterconns(icon), ...
      'numrandomization', size(statrnd,2) );

    % Step 1: compute null reference distribution
    [stat1, cfg] = clusterstat_par(cfg, statrnd, statobs);
    if numel(job.clusteralphas)*numel(job.clusterconns) == 1
      warning on
      warning('A single CDT given: non-min(p) results returned.')
      stat = stat1;
      return
    end

    if isfield(stat1,'posdistribution')
      % Step 2: convert to p-values
      [~,idx] = sort(stat1.posdistribution);
      pospvaldistribution = idx/numel(stat1.posdistribution);

      % Step 3: keep the minimal p-values only:
      posdistributionminp = min(...
        [pospvaldistribution; posdistributionminp],[],1);

      % also keep the minimal observed p-values:
      posobsminp = min([stat1.prob, posobsminp],[],2);
    end

    if isfield(stat1,'negdistribution')
      % Step 2: convert to p-values
      [~,idx] = sort(stat1.negdistribution);
      negpvaldistribution = idx/numel(stat1.negdistribution);

      % Step 3: keep the minimal p-values only:
      negdistributionminp = min(...
        [negpvaldistribution; negdistributionminp],[],1);

      % also keep the minimal observed p-values:
      negobsminp = min([stat1.prob, negobsminp],[],2);
    end

  end
end
delete(g)

% output:
cfg.clusterconns = job.clusterconns;
cfg.clusteralphas = job.clusteralphas;

stat = struct( ...
  'posdistributionminp',posdistributionminp, 'posobsminp',posobsminp, ...
  'negdistributionminp',negdistributionminp, 'negobsminp',negobsminp);

% for each voxel? can it be faster?
prb_pos = nan(nvox,1);
prb_neg = nan(nvox,1);
for ivox = 1:nvox
  prb_pos(ivox) = (sum(posdistributionminp <= posobsminp(ivox))+1) ...
    / (nrnd+1);
  prb_neg(ivox) = (sum(negdistributionminp <= negobsminp(ivox))+1) ...
    / (nrnd+1);
end

if cfg.tail==0
  % consider both tails
  % this is the probability for the most unlikely tail:
  stat.prob = min(prb_neg, prb_pos);
elseif cfg.tail==1
  % only consider the positive tail
  stat.prob = prb_pos;
elseif cfg.tail==-1
  % only consider the negative tail
  stat.prob = prb_neg;
end

end
