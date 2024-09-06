function clusterstatminp_test(Job)

%% 3D
% 1. FALSE POSITIVE CONTROL:
rng(Job.iDraw)
ArrayClass = Job.ArrayClass;
nSub = Job.nSub;
nRnd = Job.nRnd;

Job_ = [];
X = false(30,30,30);
X(2:28,2:28,2:28) = true;
Job_.inside = X;
Job_.dim = size(X);
nVox = sum(Job_.inside(:));

sqrt_nSub = sqrt(nSub);
tstat = @(x) mean(x,1) ./ std(x,[],1) .* sqrt_nSub;
Y = randn(nVox, nSub);

StatObs = tstat(Y');
StatRand = zeros(nVox, nRnd, ArrayClass); % [#voxels x #rands]

IdxToNote = fix(linspace(1, nRnd, min(nRnd,11)));
for iRnd = 1:nRnd
  Signs = repmat((2*(rand(1, nSub)<0.5)-1), [nVox 1]);
  eval(['Signs = ',ArrayClass,'(Signs);']);
  RandY = Y .* Signs;
  if gpuDeviceCount
    wait(gpuDevice)
  end
  StatRand(:,iRnd) = tstat(RandY');
  if any(iRnd==IdxToNote)
    if (iRnd>1)
      fprintf(repmat('\b', [1,LineLength]))
    end
    LineLength = fprintf('[%s:%s] Flipping signs (%2i%%)...',...
      mfilename, datestr(now,31), round(iRnd/nRnd*100));
  end
end
fprintf('\n');
Job_.clusteralphas = [0.05 0.01 0.001];
Job_.clusterconns = 6;
Job_.tail = 1;
Job_.connectivity = nan;
[stat, cfg] = clusterstatminp(Job_, StatObs, StatRand);

[~,~] = mkdir(Job.DnSim);
save(fullfile(Job.DnSim, sprintf('test_%06i.mat',Job.iDraw)), ...
  'stat','cfg','Job','StatObs')


end
