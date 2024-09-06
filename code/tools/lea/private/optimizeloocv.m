function [LambdaOpt, SSE] = optimizeloocv(X, Y, Cxx, Cxy, LambdaGrid)
% [LambdaOpt, SSE] = optimizeloocv(X, Y, Cxx, Cxy, LambdaGrid)
% Input:
%   X: {1 x nFolds} with (nSamples x nPredictors) 
%   Y: {1 x nFolds} with (nSamples x nResponses)
%   Cxx: (nPredictors x nPredictors)
%   Cxy: (nPredictors x nResponses)
%   LambdaGrid: (1 x nLambdaValues)
% Output:
%   LambdaOpt: (1 x nResponses)
%   SSE: (nLambdaValues x nResponses)

nPred = size(X{1}, 2);
nResp = size(Y{1}, 2);
nLambda = numel(LambdaGrid);
SSE = zeros(nLambda,nResp);
for iL = 1:nLambda
  L = double(LambdaGrid(iL)) * speye(nPred);
  L(1,1) = 0;
  for iSet = 1:numel(X)
    Cxx_i = Cxx - X{iSet}' * X{iSet};
    Cxy_i = Cxy - X{iSet}' * Y{iSet};
    Error = Y{iSet} - X{iSet} * ( (Cxx_i+L)\Cxy_i );
    SSE(iL,:) = SSE(iL,:) + sum( Error.^2 );
  end
end
[~, Idx] = min(SSE, [], 1);
LambdaOpt = LambdaGrid(Idx);
end
