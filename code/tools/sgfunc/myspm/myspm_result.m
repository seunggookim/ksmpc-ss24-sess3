function JOB = myspm_result(JOB)
% JOB = myspm_result(JOB)
%
% It creates SPM result reports automatically for given a directory (where you have SPM.mat)
% if there is a suprathreshold voxels/clusters,
% It generates orthogonal overlay sections for each local maximum.
%
% JOB requires:
%  .dir_glm       'Nx1' directory to save SPM results
% (.mask)         'Nx1'  filename for an explicit mask
% (.thres.desc)   'Nx1'  'FWE','none', or 'cluster'(default)
% (.thres.alpha)  [1x1]  alpha level (default=0.05)
% (.thres.extent) [1x1]  extent threshold of clusters in voxels (default=0)
% (.thres.clusterInitAlpha)  [1x1] cluster forming height threshold (default=0.001)
% (.thres.clusterInitExtent) [1x1] cluster forming extent (in voxels) threshold (default=10)
% (.fname_struct) 'Nx1' fullpath filename for background anatomical image for orthogonal slices
%                         (defulat='$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz')
% (.titlestr)     {1xNcont} Title text for SPM result report (default={'positive','negative'})
% (.dir_sum)      'Nx1' a summary directory into where you want to copy significant results
% (.append)       [1x1] whether to append results into an existing table (default=0)
% (.print)        [1x1] whether to generate figures and tables (default=1)
% (.atlas)        'Nx1' atlas to find anatomical names: 'fsl' (default) or 'spm12'
% (.MIP_only)     [1x1] prints out only MIP then exit (no slices/plots)
% (.showless)     [1x1] limit # of clusters per contrast
%
% (cc) 2015-2018, sgKIM. solleo@gmail.com, https://ggooo.wordpress.com/
if nargin < 1,  JOB=[];  end
%% Setting default parameters
%if ~isfield(JOB,'mygraph'), JOB.mygraph.y_name='y'; JOB.mygraph.x_name='x'; end
if ~isfield(JOB,'append'),  JOB.append = 0;  end
if ~isfield(JOB,'print'),   JOB.print = 1;   end
if ~isfield(JOB,'thres')||~isfield(JOB.thres,'desc'),  JOB.thres.desc='cluster';  end
if ~isfield(JOB,'dir_glm'), JOB.dir_glm=pwd; end
if ~isfield(JOB,'suppressPlot'), JOB.suppressPlot=0; end
if isfield(JOB,'masking') && ~isfield(JOB,'mask'), JOB.mask = JOB.masking; end

% set table filename
today = datestr(now,'yyyymmmdd');
fname_tab = fullfile(JOB.dir_glm,['spm_',today,'.csv']);
fmt = {   '%s',   '%s',   '%0.3f', '%s','%0.3f','%0.3f','%-0.3f','%-0.3f','%-0.3f', '%0.0f','%3.0f','%3.0f','%3.0f',  '%s','%0.2f'};
tab_fmt = cell2fmt (fmt);
if ~JOB.append && ~~JOB.print
  fid = fopen(fname_tab,'w');
  hdr_fmt = 'measure\tcontrast_name\teffect_size\tstat\tpeak\tpeakZ\tuncor_pval\tcor_pval_peak\tcor_pval_clus\tK_E\tx_mm\ty_mm\tz_mm\tpeak_strc_name\tpeak_strc_prob\n';
  fprintf(fid, hdr_fmt, JOB.thres.desc);
  fclose(fid);
end
% if ~isfield(JOB,'fname_prefix')
%   prefix='';
% else
%   prefix = [JOB.fname_prefix,'_'];
% end
if ~isfield(JOB,'model_desc')
  [~,b,~]=fileparts(JOB.dir_glm);
  JOB.model_desc=b;
end

% set figure filename
% opts=spm_print('format','pdf'); % PDF doesn't stack????
if ~isfield(JOB,'fname_spm_fig')
  JOB.fname_spm_fig = fullfile(JOB.dir_glm,[JOB.model_desc,'_spm_',today,'.ps']);
  JOB.fname_spm_fig = strrep(JOB.fname_spm_fig,'>','-gt-');
  JOB.fname_spm_fig = strrep(JOB.fname_spm_fig,'<','-lt-');
end
% 
% % delete previous figure files today (rewriting)
% if JOB.print && exist(JOB.fname_spm_fig,'file')
%   delete(JOB.fname_spm_fig);
% end

% setting a thresholding code
if strcmpi('cluster',JOB.thres.desc)
  method=3;
elseif strcmpi('FWE',JOB.thres.desc)
  method=1;
elseif strcmpi('none',JOB.thres.desc)
  method=2;
end

% setting default thresholding parameters
if ~isfield(JOB.thres,'alpha'),    JOB.thres.alpha=0.05;  end
if ~isfield(JOB.thres,'extent'),   JOB.thres.extent=0;    end
if ~isfield(JOB.thres,'override'), JOB.thres.override=0;  end
% the conventional cluster-forming alpha is double-O (0.001),
if ~isfield(JOB.thres,'clusterInitAlpha')
  JOB.thres.clusterInitAlpha  = 0.001; end
if ~isfield(JOB.thres,'clusterInitExtent')
  JOB.thres.clusterInitExtent = 10; end

% read SPM.mat to find numbers of sessions (NumSess), conditions (NumCond), and
% contrasts (NumCntrst)
load([JOB.dir_glm,'/SPM.mat']);
NumSess = size(SPM.xY.VY,2);
NumCond = size(SPM.xX.X,2)/NumSess-7;
if ~isfield(JOB,'NumCntrst')
  JOB.NumCntrst=numel(SPM.xCon);
end
NumCntrst = JOB.NumCntrst;
% set title texts
% if ~isfield(JOB,'titlestr')
for i=1:NumCntrst
  JOB.titlestr{i}=SPM.xCon(i).name;
end
% end

% find structure image for orthogonal slices
if ~isfield(JOB,'fname_struct')
  fsldir=getenv('FSLDIR');
  JOB.fname_struct = fullfile(fsldir,'data','standard','MNI152_T1_1mm.nii.gz');
end
if strcmp(JOB.fname_struct,'conmus')
  JOB.fname_struct='/home/raid2/skim/MNI152_T1_2.5mm.nii';
elseif strcmp(JOB.fname_struct,'conmus1')
  JOB.fname_struct = '/scr/vatikan4/conmus3/GLM/MNI152_T1_1mm_masked_old.nii';
elseif strcmp(JOB.fname_struct,'conmus1_masked')
  JOB.fname_struct = '/scr/vatikan4/conmus3/GLM/MNI152_T1_1mm_masked.nii';
end
[~,name1,ext1]=fileparts(JOB.fname_struct);
if strcmp(ext1,'.gz')
  if ~exist(['/tmp/',name1],'file')
    gunzip(JOB.fname_struct, '/tmp/');
  end
  JOB.fname_struct = ['/tmp/',name1];
  fprintf('Structure image file: ');
  ls(JOB.fname_struct)
end

%% Now create result reports
spm('defaults','fmri');
spm_jobman('initcfg');
JOB.minp=zeros(1,NumCntrst);
TotalNC=0;
for cntrst=1:NumCntrst %numel(JOB.titlestr) % for each contrast
  matlabbatch={};
  matlabbatch{1}.spm.stats.results.spmmat = {[JOB.dir_glm,'/SPM.mat']};
  matlabbatch{1}.spm.stats.results.conspec(1).contrasts = cntrst;
  if method == 3
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = 'none';
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = JOB.thres.clusterInitAlpha;
    matlabbatch{1}.spm.stats.results.conspec(1).extent = JOB.thres.clusterInitExtent;
  else
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = JOB.thres.desc;
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = JOB.thres.alpha;
    matlabbatch{1}.spm.stats.results.conspec(1).extent = JOB.thres.extent;
  end
  if isfield(JOB,'mask')
    matlabbatch{1}.spm.stats.results.conspec(1).mask.image.name = {[JOB.masking{1},',1']};
    matlabbatch{1}.spm.stats.results.conspec(1).mask.image.mtype = 0;
  end
  matlabbatch{1}.spm.stats.results.conspec(1).titlestr = JOB.titlestr{cntrst};
  matlabbatch{1}.spm.stats.results.units = 1;
  
  % this is what you need when you select cluster-threshold
  if method ==3
    matlabbatch{1}.spm.stats.results.print = false;
    spm_jobman('run', matlabbatch)
    
    % now read the goddamn table! :D
    TabDat = evalin('base','TabDat;');
    
    % find extent threshold for cluster-p < alpha
    Npeaks=size(TabDat.dat,1);
    K=[];
    for pi=1:Npeaks
      if (TabDat.dat{pi,3} < JOB.thres.alpha)
        K=[K TabDat.dat{pi,5}];
      end
    end
    if K
      JOB.thres.clusextent(cntrst)=min(K);
      matlabbatch{1}.spm.stats.results.conspec(1).extent = JOB.thres.clusextent(cntrst) -1;
    end
  end
  
  % result table printing
  matlabbatch{1}.spm.stats.results.print = ~~JOB.print;
  
  %% Result Reports with any given threshold
  save([JOB.dir_glm,'/result.mat'],'matlabbatch');
  spm_jobman('run', matlabbatch);
  
  %% Get all variables results report
  xSPM= evalin('base','xSPM;');
  save('xSPM.mat','xSPM') % how large?
  hReg= evalin('base','hReg;');
  TabDat = evalin('base','TabDat;');
  SPM= evalin('base','SPM;');
  
  %% Now find significant peaks
  % now read the goddamn table! :D
  
  % first just export the table (for each contrast)
  save([JOB.dir_glm,'/TabDat',num2str(cntrst),'.mat'], 'TabDat');  % for matlab
  fid = fopen([JOB.dir_glm,'/TabDat',num2str(cntrst),'.csv'],'w'); % for other applications
  fprintf(fid, cell2fmt(TabDat.hdr(1,:)));
  fprintf(fid, cell2fmt(TabDat.hdr(2,:)));
  fprintf(fid, strrep(strrep(strrep(cell2fmt(TabDat.hdr(3,:)),'\it',''),'\rm_',''),'\equiv',''));
  
  C=[7 11 3];
  Npeaks=size(TabDat.dat,1);
  COORDS={};
  PI=[];
  pvals=[];
  tmax=[];zmax=[];
  for pi=1:Npeaks
    pvals=[pvals TabDat.dat{pi,C(method)}];
    tmax=[tmax TabDat.dat{pi,9}];
    zmax=[zmax TabDat.dat{pi,10}];
    Ncols=size(TabDat.dat,2);
    for col=1:Ncols
      fprintf(fid, TabDat.fmt{col},TabDat.dat{pi,col});
      if col<Ncols
        fprintf(fid, '\t');
      else
        fprintf(fid, '\n');
      end
    end
    if (TabDat.dat{pi,C(method)} < JOB.thres.alpha)
      COORDS=[COORDS TabDat.dat{pi,end}];
      PI=[PI pi];
    end
  end
  fclose(fid);
  [~,b]=find( ~isinf(tmax),1,'first');
  cmax=tmax(b);
  JOB.cmax = cmax;
  
  % now return min P for each contrast
  JOB.minp(cntrst)=min([pvals,1]);
  NC=numel(COORDS);
  
  %% Skip the last for only MIP
  if JOB.print
    spm_print(JOB.fname_spm_fig); % add MIP + peak table in sigcluster .PDF
  end
  if isfield(JOB,'MIP_only')
    continue
  end
  
  %% Now save "significant" cluster maps (from spm code somewhere...)
  if NC
    XYZ     = xSPM.XYZ;
    Z       = round(spm_clusters(XYZ));
    num     = max(Z);
    [clusSize, ni] = sort(histc(Z,1:num), 2, 'descend');
    n       = size(ni);
    n(ni)   = 1:num;
    Z       = n(Z); % renumbering by cluster size
    
    % SANITY CHECK
    if (min(xSPM.Z)<xSPM.u) || (min(clusSize)<xSPM.k)
      error('Non-significant cluster included!')
    end
    
    fname_sigclus=['sigclus_',JOB.titlestr{cntrst}];
    fname_sigclus = strrep(fname_sigclus,'>','-gt-');
    fname_sigclus = strrep(fname_sigclus,'<','-lt-');
    fname_sigclus = strrep(fname_sigclus,' ','');
    spm_write_filtered(Z, XYZ, xSPM.DIM, xSPM.M,...
      sprintf('SPM{%c}-filtered: u = %5.3f, k = %d',xSPM.STAT,xSPM.u,xSPM.k), ...
      [fname_sigclus,'.nii']);
  end
  
  %% create orthogonal sections and reposition onto peaks
  redcmap=[ gray(64);
    [linspace(0.5667,1,14)', linspace(0,0,14)', linspace(0,0,14)'];   %  1~14 [14] dark red
    [linspace(1,1,30)', linspace(0.03333,1,30)', linspace(0,0,30)'];  % 15~44 [30] red to yellow
    [linspace(1,1,20)', linspace(1,1,20)', linspace(0.05,1,20)']      % 45~64 yellow
    ];
  
  bluecmap=[ gray(64);
    [linspace(0,0,14)', linspace(0,0,14)', linspace(0.5667,1,14)'];   %  1~14 [14] dark blue
    [linspace(0,0,30)', linspace(0.03333,1,30)', linspace(1,1,30)'];  % 15~44 [30] blue to cyan
    [linspace(0.05,1,20)', linspace(1,1,20)', linspace(1,1,20)']      % 45~64 cyan
    ];
  if ~isfield(JOB,'showlessNC') && isfield(JOB,'showless')
    JOB.showlessNC = JOB.showless;
  end
  if ~isfield(JOB,'showlessNC')
    showlessNC = 30;
  else
    showlessNC = JOB.showlessNC;
  end
  if NC>showlessNC
    warning(['Too many clusters, only ',num2str(showlessNC),' clusters are visualized']);
    NC=showlessNC;
  end
  if NC % if any cluster exists
    for ci=1:NC
      spm_sections(xSPM,hReg,JOB.fname_struct);
      spm_orthviews('reposition',COORDS{ci});
      global st
      st.vols{1}.blobs{1}.max=JOB.cmax;
      spm_orthviews('Xhirs','on');
      spm_orthviews('redraw');
      
      % odd/even index of contrast
      if mod(cntrst,2)
        colormap(redcmap);
      else
        colormap(bluecmap);
      end
      if strcmpi(xSPM.STAT,'F')
        colormap(redcmap);
      end
      
      if JOB.print
        spm_print(JOB.fname_spm_fig);
      end
      
      % and plots
      xXi = find(SPM.xCon(cntrst).c~=0);
      cfg=[]; cfg.Ic=cntrst; cfg.xXi=xXi; cfg.contrast = SPM.xCon(cntrst).c;
      cfg.origoffset=1;
      cfg.x_name = SPM.xCon(cntrst).name;
      cfg.y_name = 'y';
      if isfield(JOB,'y_name'), cfg.y_name=JOB.y_name; end
      if isfield(JOB,'markCorrThres')
        cfg.markCorrThres=JOB.markCorrThres;
        if isfield(JOB,'CorrFDRthres')
          cfg.markCorrThres.CorrFDRThres=JOB.CorrFDRthres;
        end
      end
      cfg.atlas='fsl';
      if isfield(JOB,'atlas')
        if strfind(lower(JOB.atlas),'spm')
          cfg.atlas='spm12';
        elseif strfind(lower(JOB.atlas),'fsl')
          cfg.atlas='fsl';
        end
      end
      cfg.clusterIdx = ci;
      if isfield(SPM,'Sess')
        cfg.isfmri=1;
      else
        cfg.isfmri=0;
      end
      if isfield(JOB,'suppressPlot')
        cfg.suppressPlot=JOB.suppressPlot;
      end
        
      if strcmpi(xSPM.STAT,'T')
        % SCATTER PLOTS for T-contrast
        [~,~,beta,~,STRC,thres] = myspm_graph(xSPM,SPM,hReg,cfg); % this only reads peak
        if ~isempty(thres)
          JOB.CorrFDRthres = thres;
        end
        if JOB.print
          spm_print(JOB.fname_spm_fig); % add MIP + peak table in sigcluster .PDF
        end
        %         %% fMRI: bar graph with 90% CI
        %         if any(contains({SPM.xCon.name},'Effect of interest'))
        %           cfg2 = cfg;
        %           cfg2.plotType = 'Contrast estimates and 90% C.I. effect of interest';
        %           myspm_graph(xSPM,SPM,hReg,cfg2); % this only reads peak
        %           if JOB.print
        %             spm_print(JOB.fname_spm_fig); % add MIP + peak table in sigcluster .PDF
        %           end
        %         end
        %%
      else % for F-contrast
        % only contrast estimates plot:
        cfg2 = cfg;
        cfg2.Ic = cntrst;
        cfg2.plotType = 'Contrast estimates and 90% C.I.';
        [~,~,beta,~,STRC] = myspm_graph(xSPM,SPM,hReg,cfg2); % this only reads peak
        if JOB.print
          spm_print(JOB.fname_spm_fig); % add MIP + peak table in sigcluster .PDF
        end
      end
      
      %% fMRI: PSHT?
      %         cfg2 = cfg;
      %         cfg2.plotType = 'Event-related responses';
      %         myspm_graph(xSPM,SPM,hReg,cfg2); % this only reads peak
      %         if JOB.print
      %           spm_print(JOB.fname_spm_fig); % add MIP + peak table in sigcluster .PDF
      %         end
      
      %% TODO: write correct effect size for F-contrasts
      %% TODO: write R^2 and adj. R^2
      idx_x = find(SPM.xCon(cntrst).c);
      idx_x = idx_x(1);
      % generate a summary table!
      fid = fopen(fname_tab, 'a');
      % 'y-name  x-name  effect  stattype  peak%s\tpeakZ  uncor_pval \t cor_pval(%s) \tK_E\tMNI-x_mm\tMNI-y_mm\tMNI-z_mm\tpeak-Strc-name\npeak-Strc-prob\clus-nStrc-name\nclus-Strc-prob\n'
      fprintf(fid,tab_fmt,...
        cfg.y_name, cfg.x_name, beta(idx_x), xSPM.STAT, ...
        TabDat.dat{PI(ci),9}, TabDat.dat{PI(ci),10},  ...
        TabDat.dat{PI(ci),11}, TabDat.dat{PI(ci),7}, TabDat.dat{PI(ci),3}, ...
        TabDat.dat{PI(ci),5}, ...
        TabDat.dat{PI(ci),12}(1), TabDat.dat{PI(ci),12}(2), TabDat.dat{PI(ci),12}(3), ...
        STRC.strc.name, STRC.strc.prob); %, ...
      fclose(fid);
    end % of each cluster
    TotalNC=TotalNC+NC;
  end % if any cluster exists
end % of each contrast

%% convert ps files to pdf
% myps2pdf(JOB.fname_spm_fig)

end
