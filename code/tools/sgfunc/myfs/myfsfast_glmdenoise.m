function job = myfsfast_glmdenoise(job)
%myfsfast_glmdenoise runs GLMdenoise on FS-FAST processed data
% [USAGE]
% job = myfsfast_glmdenoise(job)
%
% [INPUT]
% job (1x1) job description structure:
% - data:
%  .subid        '1xN' subject directory name
%  .subjects_dir '1xN' directory where you have the subject directory
%  .fsd          '1xN' assumes the data structure processed FS-FAST
%  .runs         [1xN] runs to smooth
%  .fwhm_mm      [1x1] full-width at half maximum in mm
%
% - design:
%  .design       '1xN' 'sound' | (...)
%  .stimdur      [1x1] duration of stim in sec
%
% - GLMdenoise options:
% (.denoisepsec) '1x5' for [model,poly,user,noise,res], 0=discard, 1=keep
%                (default: '10001' to keep only modeled responses)
% (.opt)         (1x1) a structure for GLMDENOISEDATASURF
% (.suffix_denoised) '1xN' output suffix (default: '.GLMdenoised')
%
% - output:
% (.dn_outmat)   '1xN' a parent directory where you want to have subject 
%                      directories with data & paradigms in MAT format
%                      (default: in FS-FAST structures)
%
% (cc) 2020, sgKIM, solleo@gmail.com


%% READ data:
spaces = {[job.fstemp,'.lh'],[job.fstemp,'.rh'],...
  ['mni305.',num2str(job.mni305),'mm']};
nruns = numel(job.runs);
data = cell(1,nruns);
nverts = []; ntimes = []; ind_nnz = []; tr =[]; mrihdr = {};
for irun = 1:nruns
  Y = cell(1,3);
  dn_data = fullfile(job.subjects_dir, job.subid, job.fsd, ...
    sprintf('%03i',job.runs(irun)));
  fprintf(['Reading ',dn_data,'...'])
  for s = 1:3
    mri = load_nifti(fullfile(dn_data,sprintf('%s.sm%g.%s.nii.gz',...
      job.prefix_func, job.fwhm_mm, spaces{s})));
    Y{s} = reshape(mri.vol, [], mri.dim(5))';
    mri.vol = [];
    mrihdr{irun,s} = mri;
    nverts(irun,s) = prod(mri.dim(2:4));
  end
  ntimes(irun) = mri.dim(5);
  tr(irun) = mri.pixdim(5)/1000;
  clear mri
  
  data{irun} = single([Y{1} Y{2} Y{3}]);
  isnnz = ~isnan(data{irun}) & (data{irun}~=0);
  ind_nnz = [ind_nnz; prod(isnnz)];
  clear Y
  fprintf('[#frames=%i, #voxels=%i, TR=%gs]\n', ...
    ntimes(irun), sum(nverts(irun,:)), tr(irun))
end

% Sanity check:
assert(numel(unique(sum(nverts,2))) == 1, ...
  'Data across runs are not on the same space!')
nverts = nverts(1,:);
assert(numel(unique(tr)) == 1, ...
  'TRs across runs are different!')
tr = tr(1);
job.tr = tr;
idx_nnz = find(prod(ind_nnz,1));
data = cellfun(@(x) x(:,idx_nnz)', data, 'UniformOutput', false);
fprintf('# of non-nan/zero verts = %i\n',numel(idx_nnz))


%% READ design:
design = cell(1,6);
for irun = 1:nruns
  prm = readprm([job.subjects_dir,'/',job.subid,'/',job.fsd,'/',...
    sprintf('%03i',job.runs(irun)),'/18tones.prm']);
  switch (job.design)
    case 'sound' % this is what I want to keep. right?
      % tones(1-18)
      design{irun} = zeros(ntimes(irun),1);
      times = [0:ntimes(irun)-1]*tr;
      idx_code = (prm.code >= 1) & (prm.code <= 18);
      idx_tr = dsearchn(times', prm.onset_sec(idx_code));
      design{irun}(idx_tr) = 1;
    case '18tones' % this is what I want to keep. right?
      % tones(1-18)
      design{irun} = zeros(ntimes(irun),18);
      times = [0:ntimes(irun)-1]*tr;
      for icol = 1:18
        idx_code = (prm.code == icol);
        idx_tr = dsearchn(times', prm.onset_sec(idx_code));
        design{irun}(idx_tr,icol) = 1;
      end
    case '18tones+catch+resp' % can they be teased out?
      % tones(1-18), catch(19), resp(20)
      design{irun} = zeros(ntimes(irun),20);
      times = [0:ntimes(irun)-1]*tr;
      for icol = 1:20
        idx_code = (prm.code == icol);
        idx_tr = dsearchn(times', prm.onset_sec(idx_code));
        design{irun}(idx_tr,icol) = 1;
      end
    otherwise
      error('job.design=%s NOT RECOGNIZED!', job.design)
  end
end


%% RUN GLMdenoise:
if ~isfield(job,'hrfmodel')
  job.hrfmodel = 'optimize';
end
if ~isfield(job,'hrfknobs')
  job.hrfknobs = [];
end
if ~isfield(job,'opt')
  job.opt = [];
end
if ~isfield(job,'suffix_denoised')
  job.suffix_denoised = '.GLMdenoised';
end
if ~isfield(job,'denoisespec')
  job.opt.denoisespec= '10001'; % [model,poly,exreg,pc,res]
else
  job.opt.denoisespec = job.denoisespec;
end
if ~isfield(job,'figuredir')
  job.figuredir = fullfile(job.subjects_dir, job.subid, job.fsd, ...
    ['/fig',job.suffix_denoised]);
end
surfs = myfs_readsurfs(job.fstemp, job.subjects_dir, struct('surf','inflated'));
surfs.idx_nnz = idx_nnz;
mni = [];
[mni.vol, mni.vox2ras] = ...
  load_mgh(fullfile(job.subjects_dir,'fsaverage','mri.2mm','mni305.cor.mgz'));
surfs.mri = mni;
[results,denoiseddata] = GLMdenoisedatasurf(...
  design, data, job.stimdur, job.tr, job.hrfmodel, job.hrfknobs,...
  job.opt, job.figuredir, surfs); % let's just leave the visualization out now..
clear data

% SOMEHOW THIS INCLUDES LOTS OF FALSE POSITIVES
% betas = [];
% betas.opthrf = results.modelmd{1};
% betas.prcbold = results.modelmd{2}; % in % bold change from mean
% betas.idx_nnz = idx_nnz;
% betas.nverts = nverts;
% betas.mrihdr = mrihdr(1,:);
% betas.design = design;
% betas.job = job;
% save([job.figuredir,'/betas.mat'],'betas')

results = rmfield(results,'models');
results = rmfield(results,'modelmd');
results = rmfield(results,'modelse');
save([job.figuredir,'/results.mat'],'results')

%% WRITE denoised data:
if ~isfield(job,'dn_outmat')
  for irun = 1:nruns
    dn_data = fullfile(job.subjects_dir, job.subid, job.fsd, ...
      sprintf('%03i',job.runs(irun)));
    fprintf(['Writing ',dn_data,'...'])
    
    % - put in big cell arrays (including nan/zero voxels)
    Y = nan([ntimes(irun),sum(nverts)],'single');
    if ~str2double(job.opt.denoisespec(2)) % if poly-fit was removed
      Y(:,idx_nnz) = denoiseddata{irun}' + results.meanvol';
    else
      Y(:,idx_nnz) = denoiseddata{irun}';
    end
    Y = {Y(:,1:nverts(1)), Y(:,nverts(1)+[1:nverts(2)]), ...
      Y(:,sum(nverts(1:2))+[1:nverts(3)])};
    for s = 1:3
      mri = mrihdr{irun,s};
      mri.vol = reshape(Y{s}', [mri.dim(2:4) mri.dim(5)]);
      save_nifti(mri, fullfile(dn_data,sprintf('%s.sm%g.%s%s.nii.gz',...
        job.prefix_func, job.fwhm_mm, spaces{s}, job.suffix_denoised)));
    end
    fprintf('\n')
  end
else
  dn_out = fullfile(job.dn_outmat, job.subid);
  [~,~] = mkdir(dn_out);
  for irun = 1:nruns
    % - save data
    fn_data = fullfile(dn_out, sprintf('%s_run%03i_data.mat',...
      job.subid, job.runs(irun)));
    fprintf(['Writing ',fn_data,'...'])
    data = denoiseddata{irun}';
    idx_hemi = uint8([ones(nverts(1),1); 2*ones(nverts(2),1); ...
      3*ones(nverts(3),1)]);
    hdr = struct('idx_hemi',idx_hemi, 'idx_nnz',idx_nnz, ...
      'nverts', nverts, 'ntimes', ntimes(irun), 'mrihdr', {mrihdr(irun,:)});
    save(fn_data, 'data','hdr')
  
    % - save paradigm
    prm = readprm([job.subjects_dir,'/',job.subid,'/',job.fsd,'/',...
      sprintf('%03i',job.runs(irun)),'/18tones.prm']);
    fn_prm = fullfile(dn_out, sprintf('%s_run%03i_prm.mat',...
      job.subid, job.runs(irun)));
    save(fn_prm, 'prm')
    fprintf('\n')
  end
end


%%

end


function prm = readprm(fname_prm)
fid=fopen(fname_prm);
C=textscan(fid,'%f %f %f %f %s');
prm.onset_sec = C{1};
prm.code = C{2};
prm.duration_sec = C{3};
prm.weight = C{4};
prm.txt = C{5};
end
