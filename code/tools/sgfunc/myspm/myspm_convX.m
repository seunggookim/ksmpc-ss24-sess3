function [xX,X] = myspm_convX(job)
% Creates HRf-convoluted regressors using SPM12
%
% xX = myspm_get_xX(job)
%
% job requires:
%  .tr_sec     [1x1] time of repetition in sec
%  .onset_sec  {1xR} cell array of onset times over R runs
%  .dur_sec    {1xR} cell array of durations over R runs
%  .num_frames [1xR] double array of # of frames (scans)
% (.regname)   '1xN' name of the regressor (not important; default: 'x')
% (.fMRI_T)    [1x1] number of bins within one TR (default: 16)
% (.fMRI_T0)   [1x1] order of the bin where you aligned slices (default: 8)
%
% (cc) 2020, sgKIM.

num_runs = numel(job.onset_sec);
xX = cell(1,num_runs);
if ~isfield(job,'regname'), job.regname = 'x'; end
if ~isfield(job,'fMRI_T')
  job.fMRI_T = 16;
end
if ~isfield(job,'fMRI_T0')
  job.fMRI_T0 = 8;
end
if ~iscell(job.onset_sec)
  job.onset_sec = {job.onset_sec};
end
if ~iscell(job.dur_sec)
  job.dur_sec = {job.dur_sec};
end

SPM = [];
SPM.xBF = spm_get_bf(...
  struct('T',job.fMRI_T ,'T0',job.fMRI_T0, ...
  'dt',job.tr_sec/job.fMRI_T0, 'name','hrf'));
SPM.nscan = job.num_frames;
for irun = 1:num_runs
  SPM.Sess = [];
  SPM.Sess(irun).U = struct(...
    'name',{{job.regname}},...
    'ons',job.onset_sec{irun},'dur',job.dur_sec,...
    'P',struct('name','none','h',0,'i',1));
  U = spm_get_ons(SPM, irun);
  X{irun} = spm_Volterra(U, SPM.xBF.bf, 1);
  % downsample (spm_fMRI_design:JUST DECIMATION!):
  xX{irun} = X{irun}(...
    (0:(SPM.nscan(irun) - 1))*job.fMRI_T + job.fMRI_T0 + 32,:);
  % 32 bin offset is hardcoded in spm_fMRI_design.
  %  (Perhaps something in spm_Volterra?)
end

end