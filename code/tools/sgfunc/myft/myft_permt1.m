function [stat,cfg] = myft_permt1(X1, X2, conn, k, nthreads)
% performs paired T-test (when X2 is not empty) or one-sample T-test 
% (they are actually equivalently implemented)
%
% NOTE: any high-dim data should be vectorized
%
% SYNTAX
% [stat,cfg] = myft_permt1(X1, [X2], [conn], [k], [nthreads])
%
% INPUTS
%  X1    [V x N] data in condition1 for V datapoints and N samples
% (X2)   [V x N] data in condition2
% (conn) [V x V] sparse connectivity matrix across V datapoints
% (k)    [1 x 1] # of permutation (default = 10,000)
%
% OUTPUTS
%  stat.stat [1 x V] observed stat (paired T-stat)
%  stat.randstat
%
%
% USAGE
% [stat,cfg] = myft_permt1(X0)             % without correction
%
% [stat,cfg] = myft_permt1(X1, X2)         % identical to X1-X2
%
% [stat,cfg] = myft_permt1(X1, X2, conn)   % cluster-based correction
%
% [stat,cfg] = myft_permt1(X1, X2, conn, k)  % for k times
%
% (cc) 2020, sgKIM.

[ftpath,~,~] = fileparts(which('ft_statistics_montecarlo'));
if isempty(ftpath)
  error('FieldTrip needs to be pathed')
end
if ~isfile(fullfile(ftpath,'private','clusterstat_par.m'))
  [mypath,~,~] = fileparts(mfilename('fullpath'));
  assert(~isempty(mypath),'call MYFT_PERMT1 as a function!')
  unix(['ln -s ',mypath,filesep,'cluster_stat_par.m ',...
    ftpath,filesep,'private',filesep])
end


% CHECK INPUTS:
[V,N] = size(X1);
if ~exist('X2','var') || isempty(X2)
  X2 = X1*0;
end
assert(N==size(X2,2), '# of samples are inconsistent!')
design = [
  ones(1,N) 2*ones(1,N);  % cond1 cond1 ... cond2 cond2
  repmat(1:N,[1 2]) ];    % subj1 subj2 ... subj1 subj2

% SET CONFIGURATION WITH DEFAULTS:
alpha = 0.01;
tail = 0;
  
cfg = struct('alpha',alpha, 'tail',tail, 'correcttail','prob', ...
  'ivar',1, 'uvar',2, 'statistic','depsamplesT', ...
  'numrandomization', 1e+4, 'resampling','permutation');
if exist('conn','var') && ~isempty(conn)
  assert( V==size(conn,2), '# of datapoints are inconsistent!')
  cfg.connectivity = conn;
  cfg.correctm = 'cluster';
  cfg.dimord = 'verts';
  cfg.dim = V;
end
if exist('k','var')
  cfg.numrandomization = k;
end
if V > 10
  cfg.clusterstatpar = true;
end
if exist('nthreads','var')
  cfg.nthreads = nthreads;
end

% RUN:
[stat, cfg] = ft_statistics_montecarlo(cfg, [X1 X2], design);
end
