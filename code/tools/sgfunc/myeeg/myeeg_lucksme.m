function [Stats, Bootstat, Trialstat] = myeeg_lucksme(Epoch, Job)
% Bootstrapping to estimate Luck's SME (Luck et al., 2021) of ERP components
%
% [Stats, Bootstat] = myeeg_lucksme(Epoch, Job)
%
% Epoch       <1x1 stru> EEGLAB structure, time-locked epoches
%
% Job
% .Toi_ms     <1x2 array> [min max]
% .IdxChan    <1x# array> 
% .Method     <1x# char>  'rms' (default) 'mean' 'pospeak' 'negpeak' 
% .nBoot      <1x1 array> 100000 (default)
%
% REF: Luck et al., 2021, PsychoPhysio, https://doi.org/10.1111/psyp.13793
%
% (cc) 2024, dr.seunggoo.kim@gmail.com

if not(isfield(Job,'IsSmp')) && isfield(Job,'Toi_ms')
  Job.IsSmp = (Job.Toi_ms(1) <= Epoch.times) & (Epoch.times <= Job.Toi_ms(2));
end
Job = defaultjob(struct(nBoot=100000, Method='rms'), Job, mfilename);
Data = permute(Epoch.data(Job.IdxChan, Job.IsSmp, :), [3 2 1]); % trials x time x chan (for bootstrp)
BootFunc = str2func(['erp',Job.Method]);
Bootstat = bootstrp(Job.nBoot, BootFunc, Data);
Sme = std(Bootstat, 'omitnan'); % SE is a SD of a sampling distribution
Score = feval(BootFunc, Data);
Stats = struct(Job=Job, Sme=Sme, Score=Score, SNR_dB=db(Score/Sme), ...
  BootCi=[prctile(Bootstat,2.5), prctile(Bootstat,100-2.5)]);
Trialstat = zeros(size(Data,1),1);
for j = 1:size(Data,1)
  Trialstat(j) = feval(BootFunc, Data(j,:,:));
end

end


function s = erpmean(Data)
s = mean(mean(Data,1),2);
end

function s = erprms(Data)
s = rms(mean(Data,1),2);
end

function s = erppospeak(Data)
m = mean(Data,1);
for j = size(m,3):-1:1
  pks = findpeaks(m(1,:,j), NPeaks=1);
  if isempty(pks)
    s(j) = nan; %max(m);
  else
    s(j) = pks;
  end
end
end

function s = erpnegpeak(Data)
m = mean(-Data,1);
for j = size(m,3):-1:1
  pks = findpeaks(m(1,:,j), NPeaks=1);
  if isempty(pks)
    s(j) = nan; %max(m);
  else
    s(j) = pks;
  end
end
s = -s;
end
