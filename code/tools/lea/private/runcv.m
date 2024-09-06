function Mdl = runcv(X, Y, Job)
% Mdl = runcv(X, Y, Job)
%
% Input:
%   Data: A structure array [1 x nFolds]
%   Job: A structure array [1 x 1]
%
% Output:
%   Mdl: A structure array [1 x 1] with
%     Acc: A matrix (1 x nResponses)
%     Bhat: A matrix (nFolds x nPredictors x nResponses)
tic
Data = conformdata(X, Y, Job);
[Cxx, Cxy] = findcov(Data);
nPred = size(Data(1).X, 2); 
nResp = size(Data(1).Y, 2);
Lopt = zeros(numel(Data), nResp); 
Acc = zeros(numel(Data), nResp); 
Bhat = zeros(numel(Data), nPred, nResp);

for iOuter = 1:numel(Data)
  Cxx_train = Cxx - Data(iOuter).X' * Data(iOuter).X;
  Cxy_train = Cxy - Data(iOuter).X' * Data(iOuter).Y;
  IdxTr = setdiff(1:numel(Data), iOuter);
  [Lopt(iOuter,:)] = optimizeloocv( ...
    {Data(IdxTr).X}, {Data(IdxTr).Y}, Cxx_train, Cxy_train, Job.LambdaGrid);
  [Acc(iOuter,:), Bhat(iOuter,:,:)] = evaluate( ...
    Data(iOuter).X, Data(iOuter).Y, Cxx_train, Cxy_train, Lopt(iOuter,:));
end

Mdl = struct(Lopt=Lopt, Acc=Acc, Bhat=Bhat);
logthis('DONE: took %.3f sec\n', toc)
end
