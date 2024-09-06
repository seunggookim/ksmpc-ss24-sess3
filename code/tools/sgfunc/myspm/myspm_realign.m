function JOB = myspm_realign(JOB)
% JOB = myspm_realign(JOB)
%
% unified unwarping and realignment
%
% JOB requires:
%  .fname_epi 
%
% (cc) 2015, 2019, sgKIM.  solleo@gmail.com  https://ggooo.wordpress.com

if ~isstruct(JOB)
  JOB=struct('fname_epi',JOB);
end
[p1,f1,e1] = myfileparts(JOB.fname_epi);
fname_in=[p1,'/',f1,e1];
ls(fname_in);

%% realign to MEAN IMAGE
% outputs:
% [1]  rp_${epi}.txt   : six rigid-body motion parameters [mm & rad]
% [2]  ${epi}.mat      : [4x4xT] realign transform
% [3]  ${epi}_uw.mat   : unwarping meta data
% [4]  u${epi}.nii     : unwarped/realigned image
% [5]  meanu${epi}.nii : mean image of [4]
hdr = load_untouch_header_only([p1,'/',f1,e1]);
NumFrames = hdr.dime.dim(5);
scans = [];
for t=1:NumFrames
  scans{t,1} = [p1,'/',f1,e1,',',num2str(t)];
end
matlabbatch{1}.spm.spatial.realign.estwrite.data = {scans};
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;  % because MEAN image is used in coregistration. (although RP is still w.r.t. the 1st image.. could be confusing?)
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

fname_output = [p1,'/r',f1,e1];
if ~exist(fname_output,'file')
  fname_matlabbatch=[p1,'/myspm_fmriprep4_',f1,'.mat'];
  save(fname_matlabbatch,'matlabbatch');
  spm_jobman('run', matlabbatch);
  ls(fname_output)
end
end

