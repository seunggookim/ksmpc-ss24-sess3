function [X, Y] = prepdata(Job)
% [X, Y] = prepdata(Job)
tic
assert( size(Job.FnamesX, 1) == size(Job.FnamesY, 1) )
[nSets, nFeat] = size(Job.FnamesX);
[~, nResp] = size(Job.FnamesY);
Data = {cell(nSets,1), cell(nSets,1)};
Fnames = {Job.FnamesX, Job.FnamesY};
for iSet = 1:nSets
  for iVar = 1:2
    Data{iVar}{iSet} = [];
    for iCol = 1:size(Fnames{iVar},2)
      Tbl = readtable(Fnames{iVar}{iSet, iCol});
      times = Tbl(:,1).Variables;
      vals = Tbl(:,2).Variables;

      % Align time stamps to start from zero
      fs = diff(times(1:2));
      times_zero = (0:1/fs:times(end))';
      vals_zero = interp1(times, vals, times_zero);
      vals_zero(isnan(vals_zero)) = 0;

      % apply high-pass filtering (2024a) when defined
      if Job.HighPassHz
        vals_zero = highpass(vals_zero, Job.HighPassHz, fs);
      end
      Data{iVar}{iSet} = [Data{iVar}{iSet}, vals_zero];
    end
  end
end

X = Data{1};
Y = Data{2};
logthis('DONE: took %.3f sec\n', toc)
end
