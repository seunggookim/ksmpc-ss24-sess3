function [itcWithin, itcAcross] = myft_compute_itc(timelock, toi, itc_x_conds)
% [itcWithin, itcAcross] = myft_compute_itc(timelock, [toi], [itc_x_conds])
itcWithin = [];
itcAcross = [];
conds = fieldnames(timelock);
nConds = numel(conds);
nTrl = [];
for iCond = 1:nConds
  tl = timelock.(conds{iCond});
  nTrl(iCond) = size(tl.trial,1);
end
if ~exist('toi','var')
  toi = [tl.time(1) tl.time(end)];
end
toi_idx = dsearchn(tl.time', toi')';
toi_idx = toi_idx(1):toi_idx(2);
nChans = size(tl.trial,2);
%% WITHIN EXEMPLARS
for iCond = 1:nConds
  tl = timelock.(conds{iCond});
  for iChan = 1:nChans
    trials = squeeze(tl.trial(:,iChan,toi_idx))';
    crrMtx = corr(trials);
    idx_triu = ~~(triu(ones(size(crrMtx)))-eye(size(crrMtx)));
    itcWithin(iCond,iChan) = mean(crrMtx(idx_triu));
  end
end

%% ACROSS EXEMPLARS WITHIN CONDITIONS (without randomization)
tic
for iChan = 1:nChans
  for iCond = 1:nConds
    idxCond = find(itc_x_conds == itc_x_conds(iCond));
    p = numel(idxCond);
    
    % all trials from this exemplar:
    tl = timelock.(conds{iCond});
    trials1 = squeeze(tl.trial(:,iChan,toi_idx))';
    
    % all trials from other exemplars:
    idxCondOthers = setdiff(idxCond, iCond);
    trials2 ={};
    for j = 1:numel(idxCondOthers)
      jCond = idxCondOthers(j);
      tl = timelock.(conds{jCond});
      trials2{j} = squeeze(tl.trial(:,iChan,toi_idx))';
    end
    
    % HARDCODED for 3 exemplars:
    n1 = nTrl(iCond);
    n2 = nTrl(idxCondOthers(1));
    n3 = nTrl(idxCondOthers(2));
    crrMtx1 = corr(trials1,trials2{1}); % <n1 x n2> asymmetric matrix
    crrMtx2 = corr(trials1,trials2{2}); % <n1 x n3> asymmetric matrix
    itcAcross(iCond,iChan) = mean([crrMtx1(:);crrMtx2(:)]);
  end
end
toc % 1.1097 sec

% % ACROSS CONDITIONS
% [jCond,iCond] = find(~eye(nConds));
% itcAcrossConds = [];
% for iChan = 1:nChans
%   itcX = zeros(1,numel(iCond));
%   for iPair = 1:numel(iCond)
%     tl1 = timelock.(conds{iCond(iPair)});
%     tl2 = timelock.(conds{jCond(iPair)});
%     crrMtx = corr(squeeze(tl1.trial(:,iChan,toi_idx))', squeeze(tl2.trial(:,iChan,toi_idx))');
%     itcX(iPair) = mean(crrMtx(~~(triu(ones(size(crrMtx)))-eye(size(crrMtx)))));
%   end
%   itcAcrossConds(1,iChan) = mean(itcX);
% end
% itcWithin = cat(1, itcWithin, itcAcrossConds);
% 
% %% ACROSS EXEMPLARS WITHIN CONDITIONS
% tic
% for iChan = 1:nChans
%   for iCond = 1:nConds
%     idxCond = find(itc_x_conds == itc_x_conds(iCond));
%     p = numel(idxCond);
%     
%     % all trials from this exemplar:
%     tl = timelock.(conds{iCond});
%     trials1 = squeeze(tl.trial(:,iChan,toi_idx))';
%     
%     % all trials from other exemplars:
%     idxCondOthers = setdiff(idxCond, iCond);
%     trials2 ={};
%     for j = 1:numel(idxCondOthers)
%       jCond = idxCondOthers(j);
%       tl = timelock.(conds{jCond});
%       trials2{j} = squeeze(tl.trial(:,iChan,toi_idx))';
%     end
%     
%     % HARDCODED for 3 exemplars:
%     n1 = nTrl(iCond);
%     n2 = nTrl(idxCondOthers(1));
%     n3 = nTrl(idxCondOthers(2));
%     crrMtx1 = corr(trials1,trials2{1}); % <n1 x n2> asymmetric matrix
%     crrMtx2 = corr(trials1,trials2{2}); % <n1 x n3> asymmetric matrix
%     
%     itcAcrossRand = [];
%     for iRand = 1:100
%       % randomly select k correlation coefficients from each column
%       k = round(n1*(n1-1)/2/(p-1));
%       
%       indCorr1 = false(n1*n2,1);
%       idx = randperm(numel(indCorr1));
%       indCorr1(idx(1:k)) = true;
%       
%       indCorr2 = false(n1*n3,1);
%       idx = randperm(numel(indCorr2));
%       indCorr2(idx(1:k)) = true;
%       
%       itcAcrossRand(iRand) = mean([crrMtx1(indCorr1); crrMtx2(indCorr2)]);
%     end
%     itcAcross(iCond,iChan) = mean(itcAcrossRand);
%   end
% end
% toc % 52 sec
%{ 
and this makes a REALLY small difference:
mean(abs(itc_x(:)))
ans =
    0.0079

mean(abs(itc_x_rand(:)))
ans =
    0.0080

mean(abs(itc_x(:)-itc_x_rand(:)))
ans =
   7.7973e-04
rms(itc_x(:)-itc_x_rand(:))
ans =
    0.0010
%}


end