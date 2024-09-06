function [job,cmd] = myfsl_icaaroma(job)
%Matlab wrapper for ICA-AROMA denoising package.
%
% [SYNTAX]
% [job,cmd] = myfsl_icaaroma(job)
%
% [INPUTS]
% job (1x1 structure) requires:
%  .fn_epi   '1xN char' a full path for preprocessed [normalized & smoothed] EPI
%  .fn_rp    '1xN char' a full path for rigid-motion parameters from SPM (rp_*)
% (.dn_icaaroma)  '1xN char' a full path for ICA-AROMA for a local installation (i.e., if you have multiple
% installations and if you want to run this that is not currently pathed; normally not needed) (.dn_icaaroma2p7) '1xN
% char' a full path for ICA-AROMA to run with Python 2.7
%  .runon
%
% [EXAMPLE]
% myfsl_icaaroma( struct('fn_epi','/my/path/s6wuaEPI_skip5.nii', ...
%                        'fn_rp','/my/path/rp_aEPI_skip5.nii'
%
% [OUTPUT]
% job (1x1 struc) with added fields
% cmd '1xN char'  Bash command that called Python-implemented ICA-AROMA
%
% Denoised file: for '/my/path/s6wuaEPI_skip5.nii', this wrapper will create
% '/my/path/s6wuaEPI_skip5.ica-aroma/denoised_func_data_nonaggr.nii', which is the final
% output (denoised data). Also, a symbolic link will be created for this file as
% '/my/path/s6wuaEPI_skip5_ica-aroma.nii' for convenience.
%
% (cc) 2022-11-03, dr.seunggoo.kim@gmail.com

assert(isfile(job.fn_epi), 'File "%s" not found',job.fn_epi)
assert(isfile(job.fn_rp), 'File "%s" not found',job.fn_rp)

[p1,f1,~] = myfileparts(job.fn_epi);
job.fn_mask = [p1,'/',f1,'_bet_mask.nii.gz'];
cmd = ['export FSLOUTPUTTYPE=NIFTI_GZ; bet ',p1,'/',f1,' ',p1,'/',f1,'_bet -f 0.3 -m -n -R; '];

job.dn_aroma = [p1,'/',f1,'.ica-aroma'];
if isfield(job,'dn_icaaroma') % local installation
  cmd=[cmd, ...
    ' python3 ',job.dn_icaaroma,'/ICA_AROMA.py ',...
    ' -in ',job.fn_epi,...
    ' -o ',job.dn_aroma,...
    ' -m ',job.fn_mask,...
    ' -mc ',job.fn_rp];
elseif isfield(job,'dn_icaaroma2p7') % Python 2.7
  cmd=[cmd, ...
    ' python2.7 ',job.dn_icaaroma2p7,'/ICA_AROMA.py ',...
    ' -in ',job.fn_epi,...
    ' -o ',job.dn_aroma,...
    ' -m ',job.fn_mask,...
    ' -mc ',job.fn_rp];
else % if ICA_AROMA.py is executable
  cmd=[cmd,...
    ' ICA_AROMA.py ',...
    ' -in ',job.fn_epi,...
    ' -o ',job.dn_aroma,...
    ' -m ',job.fn_mask,...
    ' -mc ',job.fn_rp];
end
if isfield(job,'dryrun') && job.dryrun
  return
end

% RUN locally or submit (and wait)
FnOut = [f1,'.ica-aroma/denoised_func_data_nonaggr.nii'];
if ~isfile([p1,'/',FnOut]) && ~isfile([p1,'/',FnOut,'.gz'])
  switch job.runon
    case 'slurm'
      slurmsh(cmd, struct(IsWait=1, Mem_GB=16)) % single precision = 8 GB
    case 'local'
      system(cmd);
    otherwise
      error('RUN THIS ON "%s"?', job.runon)
  end
end

% Ungzip only the output:
if isfile([p1,'/',FnOut,'.gz']) && ~isfile([p1,'/',FnOut])
  system(['gunzip ',p1,'/',FnOut,'.gz']);
end

% Create a relative link:
if ~isfile([p1,'/',f1,'_ica-aroma.nii'])
  system(['cd ',p1,'; ln -s ',FnOut,' ',f1,'_ica-aroma.nii']);
end

%% visualization
try
  if ~isfile([job.dn_aroma,'/classified_nonmotion_ICs.png'])
    rp = load(job.fn_rp);
    l2norm = sum([zeros(1,6); diff(rp)].^2,2);
    nvols = size(rp,1);

    % 4-D image handles
    V{1} = spm_vol(job.fn_epi);
    V{2} = spm_vol([job.dn_aroma,'/denoised_func_data_nonaggr.nii']);

    % These can be computed by reading two timepoints at once:
    tmean = cell(1,2); fdl2 = cell(1,2);
    for i = 1:2
      tmean{i} = zeros(V{i}(1).dim,'single');
      fdl2{i} = zeros(nvols,1);

      for t = 1:nvols
        imgs = spm_data_read(V{i}([t, min([t+1 nvols])]));
        fdl2{i}(t) = sum(diff(reshape(imgs,[],2)').^2);
        tmean{i} = tmean{i} + imgs(:,:,:,1);
      end
      tmean{i} = tmean{i}/nvols;
      clearvars imgs
    end

    % For these I need to read all timepoints (per slice)
    tstd = cell(1,2); rmap = cell(1,2);
    for i = 1:2
      tstd{i} = zeros(V{i}(1).dim,'single');
      rmap{i} = zeros(V{i}(1).dim,'single');
      for u3 = 1:V{i}(1).dim(3)
        Y = spm_data_read(V{i},'slice',u3);
        for u1 = 1:V{i}(1).dim(1)
          for u2 = 1:V{i}(1).dim(2)
            y = squeeze(Y(u1,u2,1,:));
            tstd{i}(u1,u2,u3) = std(y);
            rmap{i}(u1,u2,u3) = corr(l2norm, y);
          end
        end
      end
    end
  end

  % Images
  prefix = {'original','denoised'};
  for i = 1:2
    slices(tmean{i},[], struct('fname_png',[job.dn_aroma,'/',prefix{i},'_mean.png']))
    slices(tstd{i}, [], struct('fname_png',[job.dn_aroma,'/',prefix{i},'_std.png']))
    myfs_viewvol( ...
      struct('vol',tmean{i},'vox2ras',vox2ras_1to0(V{i}(1).mat)), ...
      rmap{i}, ...
      struct('xyz','axi10', 'caxis',[-.5 .5], 'thres',0.15, ...
      'fname_png',[job.dn_aroma,'/',prefix{i},'_corr-with-rp.png']))
  end

  % L2-norm plot:
  figure('position',[38   495   967   474],'visible','off')
  subplot(311)
  plot(l2norm)
  title('L2norm of frame-wise differences of rigid-body motion parameters')

  subplot(312); hold on
  plot(fdl2{1})
  plot(fdl2{2})
  legend(prefix)
  title('L2norm of frame-wise difference of images')

  subplot(313); hold on
  plot(zscore(fdl2{1}),'linewidth',2)
  plot(zscore(fdl2{2}))
  ylim0=ylim;
  ylim([-2 ylim0(2)])
  legend(prefix)
  title('L2norm of frame-wise difference of images (Z-scored)')
  export_fig([job.dn_aroma,'/plot_l2norm.png'])
  close(gcf)

catch
  warning('[%s:%s] visualization failed for "%s"', ...
    mfilename, datestr(now,31), job.fn_epi);
end

end
