function JOB = myspm_realignunwarp(JOB)
% myspm_realignunwarp(JOB)
% unwarp+realign to MEAN IMAGE
%
% JOB requires:
%  .fname_epi
% (.fname_vdm)
%
% outputs:
% [1]  rp_${epi}.txt   : six rigid-body motion parameters
% [2]  ${epi}.mat      : [4x4xT] realign transform
% [3]  ${epi}_uw.mat   : unwarping meta data
% [4]  u${epi}.nii     : unwarped/realigned image
% [5]  meanu${epi}.nii : mean image of [4]
%
% (cc) sgKIM, 2018.

if ~isfield(JOB,'overwrite'), JOB.overwrite=0; end
[p1,f1,e1]=myfileparts(JOB.fname_epi);
fname_in=[p1,'/',f1,e1];
ls(fname_in);
ru=[];
ru.data.scans=cellstr(spm_select('expand',fname_in)); % get a list of all volumes
if isfield(JOB,'fname_vdm')
  ls(JOB.fname_vdm)
  ru.data.pmscan = {[JOB.fname_vdm,',1']};
else
  ru.data.pmscan = {''};
end
ru.eoptions.quality = 0.9;
ru.eoptions.sep = 4;
ru.eoptions.fwhm = 5;
ru.eoptions.rtm = 1; % because MEAN image is used in coregistration. (although RP is still w.r.t. the 1st image.. could be confusing?)
ru.eoptions.einterp = 4;
ru.eoptions.ewrap = [0 0 0];
ru.eoptions.weight = '';
ru.uweoptions.basfcn = [12 12];
ru.uweoptions.regorder = 1;
ru.uweoptions.lambda = 100000;
ru.uweoptions.jm = 0;
ru.uweoptions.fot = [4 5];
ru.uweoptions.sot = [];
ru.uweoptions.uwfwhm = 4;
ru.uweoptions.rem = 1;
ru.uweoptions.noi = 5;
ru.uweoptions.expround = 'Average';
ru.uwroptions.uwwhich = [2 1];
ru.uwroptions.rinterp = 4;
ru.uwroptions.wrap = [0 0 0];
ru.uwroptions.mask = 1;
ru.uwroptions.prefix = 'u';
matlabbatch={};
matlabbatch{1}.spm.spatial.realignunwarp = ru;
fname_out = [p1,'/u',f1,e1];
if ~exist(fname_out,'file') || JOB.overwrite
  spm_jobman('run', matlabbatch);
  ls(fname_out);
end
end