function myspmeeg_headmodel(job)
% Prepare a headmodel using MNI mesh with fiducial points:
% "nas" [1 85 -40]
% "lpa" [-87 -12 -50]
% "rpa" [87 -12 -50]
% and save it as the 1st headmodel for the data (everything is HARD-CODED).
%
% Note: all coordinates (mm) are in the MNI space. If actual fiducial points
% stored in the MEG files are different from the anatomical points, the 
% inverse solution would be less accurate.
% 
% [USGAE]
% job = myspmeeg_headmodel(job)
%
% [INPUT]
% job structure requires:
%  .fname_meg
%
% (cc) 2020, sgKIM.

head = [];
head.D = {job.fname_meg};
head.val = 1;
head.comment = 'Anatomical fiducial points';
head.meshing.meshes.template = 1;
head.meshing.meshres = 2;

fid = [];
fid(1).fidname = 'nas';
fid(1).specification.type = [1 85 -40];
fid(2).fidname = 'lpa';
fid(2).specification.type = [-87 -12 -50];
fid(3).fidname = 'rpa';
fid(3).specification.type = [87 -12 -50];

head.coregistration.coregspecify.fiducial = fid;
head.coregistration.coregspecify.useheadshape = 0;
head.forward.eeg = 'EEG BEM';
head.forward.meg = 'Single Shell';

matlabbatch = {};
matlabbatch{1}.spm.meeg.source.headmodel = head;
spm_jobman('run',matlabbatch)


head = [];
head.D = {job.fname_meg};
head.val = 2;
head.comment = 'Fiducial coils';
head.meshing.meshes.template = 1;
head.meshing.meshres = 2;

fid = [];
fid(1).fidname = 'nas';
fid(1).specification.type = [0 85 -28];
fid(2).fidname = 'lpa';
fid(2).specification.type = [-80 -11 -47];
fid(3).fidname = 'rpa';
fid(3).specification.type = [80 -11 -47];

head.coregistration.coregspecify.fiducial = fid;
head.coregistration.coregspecify.useheadshape = 0;
head.forward.eeg = 'EEG BEM';
head.forward.meg = 'Single Shell';

matlabbatch = {};
matlabbatch{1}.spm.meeg.source.headmodel = head;
spm_jobman('run',matlabbatch)


end