function EEG = myeeg_ntprep(EEG, cfg)
% EEG = myeeg_ntprep(EEG, cfg)
%
% Preprocessing based on de Cheveigne & Arzounian (2019) Neuroimage.
% Using EEGLAB (https://sccn.ucsd.edu/eeglab/index.php) and
% NOISETOOL (http://audition.ens.fr/adc/NoiseTools/).
%
% EEG can be (1) a string of full path of .set file or (2) EEG structure
%
%
% cfg (1x1) can have the following fields:
%   .overwrite  [1x1] true=overwrite all files, false=not (default)
%
%
% (cc) 2019, sgKIM.


%% Check toolboxes needed
if ~exist('eeglab','file')
  warning(['Cannot find eeglab. You can download it from ',...
    'https://sccn.ucsd.edu/eeglab/index.php'])
  web('https://sccn.ucsd.edu/eeglab/index.php')
end
if ~exist('nt_greetings','file')
  warning(['Cannot find NoiseTool. You can download it from ',...
    'http://audition.ens.fr/adc/NoiseTools/'])
  web('http://audition.ens.fr/adc/NoiseTools/')
end


%% Set up paths for this wrapper
[mypath,~,~] = fileparts(mfilename);
addpath(genpath(mypath));
rmpath(fullfile(mypath,'private'))


%% Check EEG input
if ischar(EEG)
  [filepath,prefix,~] = fileparts(EEG);
elseif isstruct(EEG)
  filepath = EEG.filepath;
  [~,prefix,~] = fileparts(EEG.filename);
else
  error(['EEG can be (1) a string of full path of .set file or ',...
    '(2) EEGLAB structure'])
end


%% Check configuration
if ~exist('cfg','var'), cfg = []; end
if ~isfield(cfg,'overwrite'), cfg.overwrite = false; end
dn_fig = fullfile(filepath,'nt_fig');
[~,~] = mkdir(dn_fig);

%% 1. Robust detrending
if ~isfield(cfg,'nt_detrend'), cfg.nt_detrend = []; end
if ~isfield(cfg.nt_detrend,'suffix'), cfg.nt_detrend.suffix='dt'; end

% check if it's already done:
fn_out = fullfile(filepath,[prefix,'_',cfg.nt_detrend.suffix,'.set']);
if ~exist(fn_out,'file') || cfg.overwrite
  % load inputs
  EEG = pop_loadset(fullfile(filepath,[prefix,'.set']));
  
  % default parameters
  if ~isfield(cfg.nt_detrend,'order'), cfg.nt_detrend.order = 10; end
  if ~isfield(cfg.nt_detrend,'basis'), cfg.nt_detrend.basis = 'polynomials'; end
  if ~isfield(cfg.nt_detrend,'thres'), cfg.nt_detrend.thres = 3; end
  if ~isfield(cfg.nt_detrend,'niter'), cfg.nt_detrend.niter = 3; end
  if ~isfield(cfg.nt_detrend,'wsize_sec'), cfg.nt_detrend.wsize_sec = 60; end
  if ~isfield(cfg.nt_detrend,'wsize')
    cfg.nt_detrend.wsize = cfg.nt_detrend.wsize_sec*EEG.srate; % in samples
  end
  
  % run nt_detrend
  t0 = tic;
  x = EEG.data';
  fprintf([...
    '[nt_detrend] Running proc with (order=%d, basis=%s, thres=%dstd, ',...
    'niter=%d, win=%dsec)\n'], cfg.nt_detrend.order, cfg.nt_detrend.basis,...
    cfg.nt_detrend.thres, cfg.nt_detrend.niter, cfg.nt_detrend.wsize_sec)
  [x,w] = nt_detrend(x, cfg.nt_detrend.order, [], cfg.nt_detrend.basis,...
    cfg.nt_detrend.thres, cfg.nt_detrend.niter, cfg.nt_detrend.wsize);
  fprintf('\n[nt_outliers] Done in %.2f min\n',toc(t0)/60);
  
  % create figures
  fig_detrend(EEG,x,w,dn_fig,cfg,fn_out)
  
  % save outputs
  EEG.data = x';
  clear x
  EEG.nt_detrend = cfg.nt_detrend;
  EEG.nt_detrend.w = w;
  fprintf('[nt_detrend] Parameters are stored in EEG.nt_detrend\n')
  pop_saveset(EEG, fn_out);
  ls(fn_out)
end

%% 2. Outlier detection & inpainting (nt_outliers)
if ~isfield(cfg,'nt_outliers'), cfg.nt_outliers = []; end
if ~isfield(cfg.nt_outliers,'suffix'), cfg.nt_outliers.suffix='out'; end
fn_out = fullfile(filepath,[prefix,'_',cfg.nt_detrend.suffix,'_',...
  cfg.nt_outliers.suffix,'.set']);
if ~exist(fn_out,'file') || cfg.overwrite
  % load inputs
  EEG = pop_loadset(fullfile(filepath,...
    [prefix,'_',cfg.nt_detrend.suffix,'.set']));
  
  % default parameters
  if ~isfield(cfg.nt_outliers,'thres'), cfg.nt_outliers.thres = 1; end
  if ~isfield(cfg.nt_outliers,'niter'), cfg.nt_outliers.niter = 10; end
  
  % run nt_outliers
  t0 = tic;
  x = EEG.data';
  fprintf([...
    '[nt_outliers] Running proc with (thres=%d std, niter=%d)\n'], ...
    cfg.nt_outliers.thres, cfg.nt_outliers.niter)
  [w,x,h] = nt_outliers(x, EEG.nt_detrend.w, cfg.nt_outliers.thres, ...
    cfg.nt_outliers.niter);
  fprintf('[nt_outliers] Done in %.2f min\n',toc(t0)/60);
  
  % save mean difference (d) over iterations
  figure(h)
  hh=findobj(h,'type','axes');
  title(hh(2),['Filename=',EEG.filename],'interp','none')
  [~,f1,~] = fileparts(fn_out);
  export_fig(fullfile(dn_fig,[f1,'_iter.png']),'-r300')
  close(h)
  
  % create figures
  fig_outliers(EEG)
  
  % save outputs
  EEG.data = x';
  clear x
  EEG.nt_outliers = cfg.nt_outliers;
  EEG.nt_outliers.w = w;
  
  fprintf('[nt_outliers] Parameters are stored in EEG.nt_outliers\n')
  pop_saveset(EEG, fn_out);
  ls(fn_out)
  
end



%% 3. Rereferencing

%% 4. Run ICA

%% 5. Find "eye" components

%% 6. Remove line noise

%% 7. Find "occipital-alpha" components


end


function fig_detrend(EEG,x,w,dn_fig,cfg,fn_out)
%%
figure('position',[1921, 1, 946, 979], 'visible','off')
subplot(411)
imagesc(zscore(EEG.data')')
title(sprintf(...
  'Data (filename=%s, \nsrate=%d Hz, times=%d-%ds, #chans=%d, #pnts=%d)',...
  EEG.filename, EEG.srate, round(EEG.times([1 end])/1000), EEG.nbchan, EEG.pnts) ...
  ,'interp','none')
xlabel('Time [smp]'); ylabel('Chan');
cb=colorbar;  ylabel(cb, 'Z-score')

subplot(412)
imagesc(w')
colorbar;
title('Weight'); xlabel('Time [smp]'); ylabel('Chan');


subplot(413)
imagesc(zscore(EEG.data'-x)')
cb=colorbar;  ylabel(cb, 'Z-score')
title(sprintf(['Trend fit with (order=%d, basis=%s, thres=%dstd, ',...
  'niter=%d, win=%dsec)'], cfg.nt_detrend.order, cfg.nt_detrend.basis,...
  cfg.nt_detrend.thres, cfg.nt_detrend.niter, cfg.nt_detrend.wsize_sec))
xlabel('Time [smp]'); ylabel('Chan');

subplot(414)
imagesc(zscore(x)')
title('Detrended data'); xlabel('Time [smp]'); ylabel('Chan');
cb=colorbar;  ylabel(cb, 'Z-score')
[~,f1,~] = fileparts(fn_out);
export_fig(fullfile(dn_fig,[f1,'.png']),'-r300')

end

function fig_outliers()
%%
%figure('position',[1921, 1, 946, 979], 'visible','off')
figure
subplot(411)
imagesc(zscore(EEG.data')')
title(sprintf(...
  'Data (filename=%s,\n srate=%d Hz, times=%d-%ds, #chans=%d, #pnts=%d)',...
  EEG.filename, EEG.srate, round(EEG.times([1 end])/1000), EEG.nbchan, EEG.pnts) ...
  ,'interp','none')
xlabel('Time [smp]'); ylabel('Chan');
cb = colorbar;  ylabel(cb, 'Z-score')

subplot(412)
imagesc(w')
colorbar;
title('Weight'); xlabel('Time [smp]'); ylabel('Chan');

subplot(413)
imagesc(zscore(x)')
title('Inpainted data'); xlabel('Time [smp]'); ylabel('Chan');
cb = colorbar;  ylabel(cb, 'Z-score')

subplot(414)
imagesc(zscore(EEG.data'-x)')
title('Difference'); xlabel('Time [smp]'); ylabel('Chan');
cb = colorbar;  ylabel(cb, 'Z-score')


%%
[~,f1,~] = fileparts(fn_out);
export_fig(fullfile(dn_fig,[f1,'.png']),'-r300')



end