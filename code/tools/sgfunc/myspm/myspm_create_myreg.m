function JOB=myspm_create_myreg(JOB)
% creates custom regressors with convolution by HRF
% 
% JOB requires:
% .tr_sec    [1x1] TR in seconds
% .n_scans   [1x1] # of scans (volumes)
% .fname_mat '1xS' e.g. matlabfilename.mat
% .ONSETS_sec {[1 x v1], [1 x v2], ...}
% .duration_sec [1x1]
% .REGNAMES   {'1 x s1', '1 x s2', ...}
%
% (cc) 2017, sgKIM. solleo@gmail.com

% 1. create matrix
onset_tr  = round(ONSETS_sec./JOB.tr_sec);
dur_tr    = round(JOB.duration_sec/JOB.tr_sec);
K         = numel(ONSET_sec);
myregvals = zeros(n_scans,K);
for k=1:numel(K)
 idx =[];
 for j=1:dur_tr; % over the duration
  idx = [idx onset_tr(:,k)+(j-1)];
 end
 X(idx,k) = 1;
 job1=[];
 job1.TR_sec = JOB.tr_sec;
 job1.u = X(:,k);
 
 % and convolute it with HRF
 job1 = myspm_conv(job1);
 X(:,k) = job1.uc;
end

myregvals =X;
myregnames=REGNAMES;
save(fname_mat,'myregnames','myregvals');
end