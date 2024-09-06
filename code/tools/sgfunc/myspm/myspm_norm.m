function Job = myspm_norm(Job)
% Job = myspm_norm(Job)
%
% Job requires:
%  .fname_moving
%  .fname_deform (y_*)
% (.vox_mm)
% (.bbox_mm)
% (.interp) B-spline order from 0 to 6 (default=4)
%
% (cc) 2015, sgKIM, solleo@gmail.com
if nargin==0, help(mfilename); return; end
spm('Defaults','fmri');
% this is SPM12-included batch process using SPM12
a=spm('version');
if ~strcmp(a(4:5),'12'),  error(['Run ',mfilename,' on SPM12!']);  end
if ~isfield(Job,'interp'), Job.interp = 4; end % 4-spline as default

fnames={};
if iscell(Job.fname_moving)
  [path1,~,~]= fileparts(Job.fname_moving{1});
  j=1;
  for c=1:numel(Job.fname_moving)
    hdr = load_untouch_header_only(Job.fname_moving{c});
    for t=1:hdr.dime.dim(5)
      fnames{j,1} = [Job.fname_moving{c},',',num2str(t)];
      j=j+1;
    end
  end
else
  [path1,~,~]= fileparts(Job.fname_moving);
  j = 1;
  V = spm_vol(Job.fname_moving);
  for t=1:numel(V)
    fnames{j,1} = [Job.fname_moving,',',num2str(t)];
    j = j + 1;
  end
end

normalise=[];
normalise.write.subj.def = {Job.fname_deform};
normalise.write.subj.resample = fnames;
if ~isfield(Job,'bbox_mm') || isempty(Job.bbox_mm)
  normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
else
  normalise.write.woptions.bb = Job.bbox_mm;
end
if ~isfield(Job,'vox_mm')
  normalise.write.woptions.vox = hdr.dime.pixdim(2:4);
else
  normalise.write.woptions.vox = Job.vox_mm;
end
normalise.write.woptions.interp = Job.interp;

matlabbatch={};
matlabbatch{1}.spm.spatial.normalise = normalise;

save(fullfile(path1,['spm12_norm.mat']), 'matlabbatch');
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);


% quality check:
if isfield(Job, 'sanitycheck') && Job.sanitycheck
  for j = 1:numel(fnames)
    [p1,f1,e1] = myfileparts(fnames{j,1})
    fname_png = [p1,'/wua',f1,'_in_mni152.png'];
    if ~isfile(fname_png)
      fname1 = [p1,'/wua',f1,'1.nii.gz'];
      setenv('FSLOUTPUTTYPE','NIFTI_GZ');
      system(['fslroi ',p1,filesep,'wua',f1,e1,' ',...
        fname1,' 0 1']);
      fname_mni = [getenv('FSLDIR'),'/data/standard/MNI152_T1_2mm.nii.gz'];
      slices(fname1, [], ...
        struct('fname_png',fname_png,'contour',fname_mni));
    end
  end
end


end
