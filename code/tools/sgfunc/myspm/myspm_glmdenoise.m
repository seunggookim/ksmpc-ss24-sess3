function job = myspm_glmdenoise(job)
%myspm_glmdenoise a wrapper for GLMDENOISE???
%
% just READ data, CONSTRUCT design matrix from onset lists,
% SET some parameters, RUN GLMDENOISE, and WRITE denoised data.
%
% [INPUT]
% job (1x1) structure requires:
%  .fnames       {Nx1}  filenames of data of N runs
%  .stimdur_sec  [1x1]  Stim duration in sec
%  .tr_sec       [1x1]  TR in sec
%  .onset_sec    {NxP}  Onsets (sec) of P conditions over N runs
% (.figuredir)
% (.opt)
%
% (cc) 2020, sgKIM.

assert(~~exist('GLMdenoisedata','file'), 'GLMdenoisedata required')

nruns = numel(job.fnames);
nconds = size(job.onset_sec,2);
data = cell(1,nruns);
design = cell(1,nruns);
for irun = 1:nruns
  % READ data
  data{irun} = single(spm_read_vols(spm_vol(job.fnames{irun})));
  
  % CREATE design matrix
  ntimes = size(data{irun},4);
  design{irun} = zeros(ntimes, nconds);
  times_fmri = (0.5+[0:ntimes-1])'*job.tr_sec;
  for icond = 1:nconds
    idx = dsearchn(times_fmri, reshape(job.onset_sec{irun,icond},[],1));
    design{irun}(idx,:) = 1;
  end
end

if ~isfield(job,'opt')
  job.opt = struct('suffix_denoised','_GLMdenoised', ...
    'denoisespec','10001');
end
if ~isfield(job,'figuredir')
  job.figuredir = [fileparts(job.fnames{1}), job.opt.suffix_denoised];
end
if ~isfield(job,'hrfmodel')
  job.hrfmodel = 'optimize';
end


% RUN glmdenoise:
[results,denoiseddata] = GLMdenoisedata(...
  design, data, job.stimdur_sec, job.tr_sec, job.hrfmodel, [], job.opt, ...
  job.figuredir);

% SAVE fitting results:
results = rmfield(results,'models');
results = rmfield(results,'modelmd');
results = rmfield(results,'modelse');
save([job.figuredir,'/results.mat'],'results')

job.fn_out = cell(1,nruns);
for irun = 1:nruns
  % WRITE data
  [p1,f1,e1] = fileparts(job.fnames{irun});
  job.fn_out{irun} = fullfile(p1,[f1,job.opt.suffix_denoised,e1]);
  V = spm_vol(job.fnames{irun});
  for ivol = 1:size(denoiseddata{irun},4)
    Vo            = V(1);
    Vo.n          = [ivol 1];
    Vo.fname      = job.fn_out{irun};
    Vo.descrip    = sprintf('GLMdenoised');
    if ~str2double(job.opt.denoisespec(2)) % if poly-fit was removed
      interp = results.meanvol;
    else
      interp = 0;
    end
    img = denoiseddata{irun}(:,:,:,ivol) + interp;
    eval(['img = ',spm_type(4),'(img);']); % put the precision back.
    spm_write_vol(Vo, img);
  end
end

end