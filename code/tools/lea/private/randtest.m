function Rnd = randtest(X, Y, Mdl, Job)
% Rnd = randtest(Mdl, Job)

tic
AccRnd = cell(Job.nRands, 1);
for iRnd = 1:Job.nRands % for parfor?
  AccRnd{iRnd,1} = randthis(X, Y, Mdl, Job);
end

AccRnd = cell2mat(AccRnd);
PvalUnc = mean([mean(Mdl.Acc,1); AccRnd] >= mean(Mdl.Acc,1));
[~,~,PvalFdr] = fdr(PvalUnc);

Rnd = struct(PvalUnc=PvalUnc, PvalFdr=PvalFdr, AccRnd=AccRnd);
logthis('DONE: took %.3f sec\n', toc)

end

function Acc = randthis(X, Y, Mdl, Job)
RandX = X;
for i = 1:numel(X)
  RandX{i} = randomize_phase(X{i});
end
Data = conformdata(RandX, Y, Job);
[Cxx, Cxy] = findcov(Data);
nResp = size(Data(1).Y, 2);
Acc = zeros(numel(Data), nResp); 
% Bhat = zeros(numel(Data), nReg, nResp);

for iOuter = 1:numel(Data)
  Cxx_train = Cxx - Data(iOuter).X' * Data(iOuter).X;
  Cxy_train = Cxy - Data(iOuter).X' * Data(iOuter).Y;
  [Acc(iOuter,:)] = evaluate( ...
    Data(iOuter).X, Data(iOuter).Y, Cxx_train, Cxy_train, Mdl.Lopt(iOuter,:));
end

Acc = mean(Acc,1);
end
