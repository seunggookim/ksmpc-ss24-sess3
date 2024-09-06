function myspm_prepare_vdm(shortmag, phase, TEs_msec, trgEPI1, TotalReadoutTime_msec, fname_t1w)
% myspm_prepare_vdm(shortmag, phase, TEs, trgEPI, TotalReadoutTime, fname_t1w)
%
% shortmag: short TE magnitude
% phase: phase difference
% TEs : [te1 te2]
% TotalReadoutTime [msec]: Total readout time of the >>TARGET EPI<< (not the fieldmap!)
% trgEPI1: epi to unwarp >>>NOTE<<<: will create a warped version but only
% with the first volume of the 4-D image. So please feed in just the first
% volume (or mean volume) so that it won't create a file named
% "ua${EPI}.nii" but with only the first volume of the 4-D timeseries by
% mistake.
%
% (cc) 2014, 2019, sgKIM, solleo@gmail.com

disp(['> TotalReadoutTime : ',num2str(TotalReadoutTime_msec),' msec']);
disp(['> ShortEchoTime    : ',num2str(TEs_msec(1)),' msec']);
disp(['> LongEchoTime     : ',num2str(TEs_msec(2)),' msec']);

%% setup matlabbatch
[path1,f1,e1] = myfileparts(phase);
phase=[f1,e1];
[~,f1,e1] = myfileparts(shortmag);
shortmag=[f1,e1];

subj=[];
subj.phase{1}     = [fullfile(path1,phase),',1'];
subj.magnitude{1} = [fullfile(path1,shortmag),',1'];
subj.defaults.defaultsval.et = TEs_msec;
subj.defaults.defaultsval.maskbrain = +1;
subj.defaults.defaultsval.blipdir   = -1;
subj.defaults.defaultsval.tert      = TotalReadoutTime_msec;
subj.defaults.defaultsval.epifm     =  1;
subj.defaults.defaultsval.ajm       =  0;

subj.defaults.defaultsval.uflags.method = 'Mark3D';
subj.defaults.defaultsval.uflags.fwhm   = 10;
subj.defaults.defaultsval.uflags.pad    =  0;
subj.defaults.defaultsval.uflags.ws     =  1;

subj.defaults.defaultsval.mflags.template{1} = fullfile(fileparts(which('spm')),'templates','T1.nii');
if isempty(dir(subj.defaults.defaultsval.mflags.template{1}))
 subj.defaults.defaultsval.mflags.template{1} = fullfile(fileparts(which('spm')),'toolbox','OldNorm','T1.nii');
end
subj.defaults.defaultsval.mflags.fwhm     = 5;
subj.defaults.defaultsval.mflags.nerode   = 2;
subj.defaults.defaultsval.mflags.ndilate  = 4;
subj.defaults.defaultsval.mflags.thresh   = 0.5;
subj.defaults.defaultsval.mflags.reg      = 0.2;

[path1,f1,e1]=fileparts(trgEPI1);
if isempty(path1), path1=pwd; end
trgEPI1=[f1,e1];

subj.session.epi{1} = [fullfile(path1,trgEPI1),',1'];
subj.matchvdm       = 1;
subj.sessname       = 'session';
subj.writeunwarped  = 1;
subj.anat           = '';
subj.matchanat      = 0;
if exist('fname_t1w','var')
 if ~isempty(fname_t1w)
  [p2,f2,e2] = myfileparts(fname_t1w);
  fn = [tempname,e2];
  copyfile(fname_t1w, fn);
  subj.anat = {[fn,',1']};
  subj.matchanat    = 1;
 end
end

matlabbatch={};
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj = subj;

spm_jobman('initcfg')
spm_jobman('run', matlabbatch);
end
