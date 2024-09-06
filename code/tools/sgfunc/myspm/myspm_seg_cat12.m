function job = myspm_seg_cat12(job, runmode)
% job = myspm_seg_cat12(job, [runmode]) runs CAT12 segmentation, which uses
% optimized Shooting normalization & a bit of custom approches.
%
% job requires:
%  .fname_t1w  [string]
% (.iseastern) [1x1] default=0 (if 1 use East-Asian template for affine transform)
% (.mw)        [1x5] write (true/false) normalized & modulated images for
%                    GM/WM/CSF/HWM/SL tissue types
% (.native)
% (.norm)
% (.surface)   [1x1] 0=Nope, 1=standard, 5=low_res_only_for_preview?
%
% if isMp2rage=1, then the function will use __original image__ as the input instead of intensity-bias
% corrected one (i.e., m*)
%
% (cc) 2019, sgKIM, solleo@gmail.com

if nargin == 0, help(mfilename); return; end
if ~isfolder(fullfile(spm('dir'),'toolbox','cat12'))
  error('This function requires CAT12 Toolbox installed. Download it from http://www.neuro.uni-jena.de/cat/')
end

if ~isstruct(job) && ischar(job), job=struct('fname_t1w',job); end
if ~isfield(job,'iseastern') || (isfield(job,'iseastern') && ~job.iseastern)
  AFFREGTRG = 'mni';
else
  AFFREGTRG = 'eastern';
end
if ~isfield(job,'surface'), job.surface = false; end
if ~isfield(job,'norm'), job.norm=0; end
if ~isfield(job,runmode)
  switch job.runmode
    case {'ss'}
      job.mw = zeros(1,5);
      job.native = [1 1 1 0 0];
      job.norm = 1;
  end
else
  runmode = '';
end

% check input file:
ls(job.fname_t1w);
% [p2,f2,e2] = myfileparts(job.fname_t1w);
% dn_cat = fullfile(p2,'/cat');
% [~, ~] = mkdir(dn_cat);
% % cd(dn_cat)
% try
%   unix(['ln -sf ',fullfile(p2,[f2,e2]),' ',dn_cat])
% catch ME
%   copyfile(job.fname_t1w, dn_cat)
% end

% run cat12:
matlabbatch = {};
matlabbatch{1}.spm.tools.cat.estwrite.data = {job.fname_t1w};
matlabbatch{1}.spm.tools.cat.estwrite.nproc = 0;
matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = {fullfile(spm('dir'),'tpm','TPM.nii')};
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = AFFREGTRG;
matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.opts.accstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.APP = 1070;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.NCstr = -Inf;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.LASstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.gcutstr = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.cleanupstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.WMHC = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.SLC = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.restypes.fixed = [1 0.1];
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.shootingtpm = {'/home/sgk/Documents/MATLAB/spm12/toolbox/cat12/templates_1.50mm/Template_0_IXI555_MNI152_GS.nii'};
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.regstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox = 1.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtres = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.scale_cortex = 0.7;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.add_parahipp = 0.1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.close_parahipp = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.ignoreErrors = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.verb = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.print = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.surface = false;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.aal = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.mori = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.anatomy = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = job.native(1);
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = job.mw(1);
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = job.native(2);
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = job.mw(2);
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native = job.native(3);
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod = job.mw(3);
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.native = job.native(4);
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.mod = job.mw(4);
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.native = job.native(5);
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.mod = job.mw(5);
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobian.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [1 0];
matlabbatch{1}.spm.tools.cat.estwrite.output.surface = job.surface;

% save test.mat matlabbatch

% % make skull-stripped image:
% fname_out = [p2,filesep,'bm',f2,e2];
% if ~isfile(fname_out)
  %   cat12('expert')
  cat_defaults
  global cat
  cat.extopts.expertgui    = 1;     % control of user GUI:   0 - common user modus with simple GUI; 1 - expert modus with extended GUI; 2 - developer modus with full GUI
  cat.extopts.subfolders   = 0;     % use subfolders such as mri, surf, report, and label to organize your data
  cat.extopts.print        = 0;     % display and print out pdf-file of results: 0 - off, 2 - volume only, 2 - volume and surface (default)
  spm_jobman('initcfg');
  spm_jobman('run',matlabbatch);
    spm('quit')
%   
%   Y = cell(1,3);
%   for c = 1:3
%     V = spm_vol_nifti([dn_cat,filesep,'p',num2str(c),f2,e2]);
%     if strcmp(runmode, 'ss')
%       copyfile([dn_cat,filesep,'p',num2str(c),f2,e2],...
%         [p2,filesep,'c',num2str(c),f2,e2])
%     end
%     [Y{c},~] = spm_read_vols(V);
%   end
%   V = spm_vol_nifti([dn_cat,filesep,'m',f2,e2]);
%   
%   manat_uncut=spm_read_vols(V);
%   manat = manat_uncut;
%   manat(~(Y{1}|Y{2}|Y{3}))=0;
%   V.fname=fname_out;
%   V.descrip='Skull stripped';
%   spm_write_vol(V,manat);
%   
%   try
%     % quality check:
%     fname_png=[p2,filesep,'bm',f2,'.png'];
%     slices(manat_uncut, struct('fname_png',fname_png), double(Y{1}|Y{2}|Y{3}))
%   end
% end
% 
% if job.norm
%   [p3,f3,e3] = myfileparts(fname_out);
%   fname_out = [p3,filesep,'w',f3,e3];
%   if ~exist(fname_out,'file')
%     copyfile([dn_cat,filesep,'wm',f2,e2], fname_out)
%     try
%       % quality check:
%       fname_png=[p3,filesep,'w',f3,'_in_mni152.png'];
%       fname_mni=fullfile(getenv('FSLDIR'),'data','standard','MNI152_T1_1mm_brain.nii.gz');
%       slices(fname_out, fname_mni, struct('fname_png',fname_png))
%     end
%   end
% end
% 
% if strcmp(runmode, 'ss')
%   copyfile([dn_cat,filesep,'y_',f2,e2],...
%     [p2,filesep])
% end

end
