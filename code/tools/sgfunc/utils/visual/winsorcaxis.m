function [caxis, caxis_unsigned, class] = winsorcaxis(numvals, doTEST)
% [caxis, caxis_unsigned] = windorcais(numvals)
%

if exist('doTEST','var') && doTEST
  TEST
  return
end

%{
positive continuous values: Winsorized range
signed continous values: mirrored winsorized range
binray: [0 1]
class labels (integers): [min max] (why not -1 and +1?)
single value: [thatValue-1 thatValue]
%}

numvals = numvals(:);
numvals(isnan(numvals)) = [];
numvals(isinf(numvals)) = [];
minval = min(numvals);
maxval = max(numvals);
p1 = prctile(numvals,1);
p99 = prctile(numvals,99);
nlevels = numel(unique(numvals));

% WHAT is LABEL? 
% 1. all values are integers (positive or negative)
% 2. and..?
isLabel = isequal(unique(round(numvals)),unique(numvals));

if (nlevels==0)
  caxis = [0 1];
  caxis_unsigned = [0 1];
  class = 'null';
elseif (nlevels==1)  % if single-valued
  caxis = [minval-1, minval];
  caxis_unsigned = [minval-1, minval];
  class = 'single-value';
elseif isLabel
  caxis = [minval, maxval];
  caxis_unsigned = [minval, maxval];
  class = 'labels';
elseif (minval*maxval)>=0 % same signed
  caxis = [p1, p99];
  caxis_unsigned = [p1, p99];
  class = 'same-signed-continuous-values';
else % signed
  caxis = [min(-abs(p1),-abs(p99)) max(abs(p1),abs(p99))];
  caxis_unsigned = [p1, p99];
  class = 'signed-continuous-values';
end

end


%% UNIT TESTING
function TEST
rng('default')
x = round(randn(100,1)*1000);
y = 0.1:0.1:1;

testnames = {'empty','nan','zeros','ones','singleval','binary','posint'};
tests     = {[],   nan,  [0 0 0],[1 1 1],[.5 .5 .5],[0 1 1],[1 2 3]};
targets_a = {[0 1],[0 1],[-1 0], [0 1],  [-.5 .5],  [0 1],  [1 3]};
targets_b = {[0 1],[0 1],[-1 0], [0 1],  [-.5 .5],  [0 1],  [1 3]};

testnames = [testnames , {'signedint','poscont','negcont','signedcont'}];
tests     = [tests,     {[x]         ,[y]      ,[-y]     ,[-1:.1:1]}];
targets_a = [targets_a, {[-2944 3578],[.1 1]   ,[-1 -.1] ,[-1 1]}];
targets_b = [targets_b, {[-2944 3578],[.1 1]   ,[-1 -.1] ,[-1 1]}];

for itest = 1:numel(tests)
  [a,b,c] = winsorcaxis(tests{itest});
  disp(['for ',testnames{itest},':', c])
  assert(isequal(a, targets_a{itest}))
  assert(isequal(b, targets_b{itest}))
end


end