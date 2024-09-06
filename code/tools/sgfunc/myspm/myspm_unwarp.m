function JOB = myspm_unwarp(JOB)
% JOB = myspm_unwarp(JOB)
%
% unified unwarping and realignment
%
% JOB requires:
%  .fname_epi 
% (.fname_mag, fname_pha) with json files!
%  .fname_t1w  for VDM preparation
%
% (cc) 2015, 2019, sgKIM.  solleo@gmail.com  https://ggooo.wordpress.com

[p1,f1,e1] = myfileparts(JOB.fname_epi);

%% preparing VDM
if isfield(JOB,'fname_mag') && isfield(JOB,'fname_pha')
  if ~isfield(JOB,'TEs_fmap')
    [p5,f5,~] = myfileparts(JOB.fname_mag);
    try
      json1 = jsondecode(fileread([p5,'/',f5,'.json']));
      [p5,f5,~] = myfileparts(JOB.fname_pha);
      json2 = jsondecode(fileread([p5,'/',f5,'.json']));
      if ~isfield(json2,'EchoTime') && isfield(json2,'EchoTime2')
        JOB.TEs_fmap = [json1.EchoTime json2.EchoTime2]*1000;
      else
        JOB.TEs_fmap = [json1.EchoTime json2.EchoTime]*1000;
      end
      %         JOB.totalreadout_msec = 1/json2.PixelBandwidth*1000;
    catch
      error(['Cannot find json files for fieldmap. ',...
        'Enter .TEs_fmap and .totalreadout_msec manually']);
    end
  end
  if ~isfield(JOB,'totalreadout_msec')
    if strcmp(f1(1),'a')
      json = jsondecode(fileread(fullfile(p1,[f1(2:end),'.json'])));
    else
      json = jsondecode(fileread(fullfile(p1,[f1,'.json'])));
    end
    JOB.totalreadout_msec = json.TotalReadoutTime * 1000; % THIS IS of the EPI to unwarp!
  end
  
  % Create EPI-specific instance:
  [p0,f0,e0]= fileparts(JOB.fname_mag);
  fn = fullfile(p0,[f0,'_',f1,e0]);
  system(['ln -s ',JOB.fname_mag,' ',fn]);
  JOB.fname_mag = fn;
  
  [p0,f0,e0]= fileparts(JOB.fname_pha);
  fn = fullfile(p0,[f0,'_',f1,e0]);
  system(['ln -s ',JOB.fname_pha,' ',fn]);
  JOB.fname_pha = fn;
  
  [~,f_ph,e_ph] = myfileparts(JOB.fname_pha);
  JOB.fname_vdm = [p0,'/vdm5_sc',f_ph,e_ph];
  if ~isfile(JOB.fname_vdm)
    % CROP the first volume to fit the fieldmap
    V = spm_vol(JOB.fname_epi);
    V1 = V(1);
    V1.fname = [p1,'/',f1,'1',e1];
    spm_write_vol(V1, spm_read_vols(V(1)));
    
    % CREATE voxel-displacement map (to correct B0-related distortion)
    myspm_prepare_vdm(JOB.fname_mag, JOB.fname_pha, JOB.TEs_fmap, ...
      [p1,'/',f1,'1',e1], JOB.totalreadout_msec, JOB.fname_t1w);
    
    % visualize unwarping results:
    [p2,f2,e2] = myfileparts(JOB.fname_t1w);
    if ~isfile([p2,'/om',f2,e2])
      myspm_coreg(struct('fname_moving',[p2,'/m',f2,e2], ...
        'fname_fixed',[p1,'/',f1,'1',e1]))
    end
    compare_unwarped([p1,'/',f1,'1',e1], [p1,'/u',f1,'1',e1], ...
      [p2,'/om',f2,e2])
  end
  ls(JOB.fname_vdm)
end
%% unwarp+realign to MEAN IMAGE
% outputs:
% [1]  rp_${epi}.txt   : six rigid-body motion parameters [mm & rad]
% [2]  ${epi}.mat      : [4x4xT] realign transform
% [3]  ${epi}_uw.mat   : unwarping meta data
% [4]  u${epi}.nii     : unwarped/realigned image
% [5]  meanu${epi}.nii : mean image of [4]
V = spm_vol([p1,'/',f1,e1]);
NumFrames = numel(V);
realignunwarp1=[];
for t=1:NumFrames
  realignunwarp1.data.scans{t,1} = [p1,'/',f1,e1,',',num2str(t)];
end
if isfield(JOB,'fname_vdm')
  ls(JOB.fname_vdm)
  realignunwarp1.data.pmscan = {[JOB.fname_vdm,',1']};
else
  realignunwarp1.data.pmscan = {''};
end
realignunwarp1.eoptions.quality = 1;
realignunwarp1.eoptions.sep = 4;
if isfield(JOB,'realign_sep')
  realignunwarp1.eoptions.sep = JOB.realign_sep;
end
realignunwarp1.eoptions.fwhm = 5;
realignunwarp1.eoptions.rtm = 1; % because MEAN image is used in coregistration. (although RP is still w.r.t. the 1st image.. could be confusing?)
realignunwarp1.eoptions.einterp = 2;
realignunwarp1.eoptions.ewrap = [0 0 0];
realignunwarp1.eoptions.weight = '';
realignunwarp1.uweoptions.basfcn = [12 12];
realignunwarp1.uweoptions.regorder = 1;
realignunwarp1.uweoptions.lambda = 100000;
realignunwarp1.uweoptions.jm = 0;
realignunwarp1.uweoptions.fot = [4 5];
realignunwarp1.uweoptions.sot = [];
realignunwarp1.uweoptions.uwfwhm = 4;
realignunwarp1.uweoptions.rem = 1;
realignunwarp1.uweoptions.noi = 5;
realignunwarp1.uweoptions.expround = 'Average';
realignunwarp1.uwroptions.uwwhich = [2 1];
realignunwarp1.uwroptions.rinterp = 4;
realignunwarp1.uwroptions.wrap = [0 0 0];
realignunwarp1.uwroptions.mask = 1;
realignunwarp1.uwroptions.prefix = 'u';
matlabbatch={};
matlabbatch{1}.spm.spatial.realignunwarp = realignunwarp1;
fname_in = [p1,'/',f1,e1];
ls(fname_in);
fname_output = [p1,'/meanu',f1,e1];
if ~exist(fname_output,'file') || JOB.overwrite
  fname_matlabbatch=[p1,'/myspm_fmriprep4_',f1,'.mat'];
  save(fname_matlabbatch,'matlabbatch');
  spm_jobman('run', matlabbatch);
  ls(fname_output)
end
end

