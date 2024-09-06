function job = myspm_viewsurf_glms(job)
% job = myspm_viewsurf_glms(job)
% job requires:
%  .dir_glm  : path for SPM results directory
%  .cfg      : cfg structure to pass to MYFS_VIEW
% (cc) 2019, 2021, sgKIM.

NOFIGURE = 1; % For normal usage; 0= for debugging

cd(job.dir_glm)
load('SPM.mat','SPM')
if ~isfield(job,'fn_bigfig')
  job.fn_bigfig = [pwd,'/bigfig.png'];
end
if ~isfield(job,'cfg'), job.cfg = []; end
%% (FIND mask if not given)
if ~isfield(job.cfg, 'masks')
  [~,~,job.cfg.masks] = myfs_mni2fs('mask.nii', struct('fsavg','fsaverage6',...
    'nofigure',1, 'interp','nearest'));
end

%%
fname_png = {};
%% Mean Image
if isfield(job,'mean')
  % FIND constant term:
  if ~isfield(SPM,'Sess')
    meanimage = 'beta_0001';
  else
    error('How can I find it?')
    meanimage = '???';
  end
  [~, surfs, Y] ...
    = fsss_mni_to_fsavg([meanimage,'.nii'], ...
    'fsaverage6', struct('nofigure',NOFIGURE,'interpopt','nearest',...
    'mapmethod','nnf'));
  
  
  cfg = job.cfg;
  cfg.fname_png = [tempname,'.png'];
  fname_png = [fname_png cfg.fname_png];
  cfg.colorbartitle = 'Mean';
  
  myfs_view(surfs, Y, cfg)
end
%%
if ~isfield(job,'idxCntrst')
  job.idxCntrst = 1:2:numel(SPM.xCon);
end
for i = job.idxCntrst
  fn_src = sprintf('%s_%04i.nii',job.prefix,i);
  if ~exist(fn_src,'file')
    warning('%s not found. skipping it.',fn_src)
    continue
  end
  [~, surfs, Y] = fsss_mni_to_fsavg(fn_src, ...
    'fsaverage6', struct('nofigure',NOFIGURE,'interpopt','nearest',...
    'mapmethod','nnf'));
  
  cfg = job.cfg;
  if ~isfield(cfg,'fname_png')
    cfg.fname_png = [tempname,'.png'];
  end
  fname_png = [fname_png cfg.fname_png];
  if ~isfield(cfg,'colorbartitle')
    cfg.colorbartitle = ['Effect of ',SPM.xCon(i).name];
  end
  if isfield(job,'caxis')
    cfg.caxis = job.caxis;
  end
  if ~isfield(job,'thresp')
    job.thresp = 0.05;
  end
  if strcmp(SPM.xCon(i).STAT,'T')
    df = SPM.xX.erdf;
    cfg.colorbarxlabel = ['T(',num2str(round(df)),') | unc-p<',...
      num2str(job.thresp)];
    cfg.thres = icdf('t',1-job.thresp,df);
    if job.thresp == 1
      cfg.thres = 0;
    end
    if ~isfield(job,'caxis')
      cfg.caxis = [-10 10];
    else
      cfg.caxis = job.caxis;
    end
  else
    df = [SPM.xCon(i).eidf SPM.xX.erdf];
    cfg.thres = icdf('f', 1-job.thresp, df(1), df(2));
    if job.thresp == 1
      cfg.thres = 0;
    end
    cfg.colorbarxlabel = ['F(',num2str(round(df(1))),...
      ',',num2str(round(df(2))),' | unc-p<',num2str(job.thresp)];
    if ~isfield(job,'caxis')
      cfg.caxis = [0 20];
    else
      cfg.caxis = job.caxis;
    end
  end
  
  if isfield(job,'colorbarxlabelsuffix')
    cfg.colorbarxlabel = [cfg.colorbarxlabel job.colorbarxlabelsuffix];
  end
  
  Y{1}(Y{1}==0) = nan; Y{2}(Y{2}==0) = nan;
  if isfield(job,'masks')
    cfg.masks = job.masks';
  end
  
  if isfield(job,'sigclus') && job.sigclus
    if strcmp(job.prefix(end),'T')
      sign = {'+','-'};
      sigclus = cell(1,2);
      for jSign = 1:2
        fn_sigclus = ['sigclus_',sign{jSign},SPM.xCon(i).name(2:end),'.nii'];
        if exist(fn_sigclus,'file')
          [~, surfs, sigclus{jSign}] = fsss_mni_to_fsavg(fn_sigclus, ...
            'fsaverage6',struct('nofigure',NOFIGURE,'interpopt','nearest',...
            'mapmethod','nnf'));
        else
          warning('%s not found',fn_sigclus)
          sigclus{jSign} = {false(size(Y{1})), false(size(Y{2}))};
        end
      end
      cfg.masks = {};
      cfg.masks{1} = ~~(sigclus{1}{1} + sigclus{2}{1});
      cfg.masks{2} = ~~(sigclus{1}{2} + sigclus{2}{2});
      
    elseif strcmp(job.prefix(end),'F')
      fn_sigclus = ['sigclus_',SPM.xCon(i).name(1:end),'.nii'];
      [~, surfs, sigclus] = fsss_mni_to_fsavg(fn_sigclus, ...
        'fsaverage6',struct('nofigure',NOFIGURE,'interpopt','nearest',...
        'mapmethod','nnf'));
      cfg.masks = {~~sigclus{1}, ~~sigclus{2}};
      
    else
      error('job.prefix=%s NOT FOR SIGNIFICANT CLUSTERS!',job.prefix)
    end
  end
  
  fsss_view(surfs, Y, cfg)
end
imageconcat(fname_png, job.fn_bigfig, 1, struct('scale',0.25))

end