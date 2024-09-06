function [Acc, Bhat] = evaluate(Xi, Yi, Cxx, Cxy, LambdaOpt)
% [Acc, Bhat] = evaluate(Xi, Yi, Cxx, Cxy, LambdaOpt)
% Input:
%   Xi: A matrix (nSamples x nRegressors) 
%   Yi: A matrix (nSamples x nResponses)
%   Cxx: A matrix  (nPredictors x nPredictors)
%   Cxy: A matrix (nPredictors x nResponses)
%   LambdaOpt: A vector (1 x nResponses)
% Output:
%   Acc: A vector (1 x nResponses)
%   Bhat: A matrix (nPredictors x nResponses)

UniqueLambdas = unique(LambdaOpt);
[~, nPred] = size(Xi);
[~, nResp] = size(Yi);
Bhat = nan(nPred, nResp, 'double');
for iL = 1:numel(UniqueLambdas)
  IdxResp = LambdaOpt == UniqueLambdas(iL);
  L = double(UniqueLambdas(iL)) * speye(nPred);
  L(1,1) = 0;
  Bhat(:,IdxResp) = (Cxx+L)\Cxy(:,IdxResp);
end
Acc = corrvec(Yi, Xi*Bhat);
end
