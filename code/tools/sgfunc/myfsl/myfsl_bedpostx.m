function JOB = myfsl_bedpostx(JOB)
% myfsl_probtrackx2(JOB)
%
% JOB requires:
% .fn_dwi
% .fn_nodifmask

[p1,f1,e1] = myfileparts(JOB.fn_dwi);
unix(['ln -sf ',JOB.fn_dwi,' data',e1]);
unix(['ln -sf ',p1,'/',f1,'.bval ',p1,'/bvals']);
unix(['ln -sf ',p1,'/',f1,'.bvec ',p1,'/bvecs']);
[~,~,e2] = myfileparts(JOB.fn_nodifmask);
unix(['ln -sf ',JOB.fn_nodifmask,' ',p1,'/nodif_brain_mask',e2]);

%% BEDPOSTX
%{
$ bedpostx

Usage: bedpostx <subject directory> [options]

expects to find bvals and bvecs in subject directory
expects to find data and nodif_brain_mask in subject directory
expects to find grad_dev in subject directory, if -g is set
options (old syntax)
-n (number of fibres per voxel, default 3)
-w (ARD weight, more weight means less secondary fibres per voxel, default 1)
-b (burnin period, default 1000)
-j (number of jumps, default 1250)
-s (sample every, default 25)
-model (Deconvolution model. 1: with sticks, 2: with sticks with a range of diffusivities (default), 3: with zeppelins)
-g (consider gradient nonlinearities, default off)


ALTERNATIVELY: you can pass on xfibres options onto directly bedpostx
 For example:  bedpostx <subject directory> --noard --cnonlinear
 Type 'xfibres --help' for a list of available options 
 Default options will be bedpostx default (see above), and not xfibres default.

Note: Use EITHER old OR new syntax.
%}
setenv('FSLOUTPUTTYPE','NIFTI_GZ');
if ~exist([p1,'/merged_f1samples.nii.gz'],'file')
 unix(['bedpostx ',p1,' -n 3 -w 1 -b 1000 -j 1250 -s 25 -model 2 '])
end