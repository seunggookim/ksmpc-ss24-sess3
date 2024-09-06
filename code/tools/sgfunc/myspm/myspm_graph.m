function [Y,y,beta,Bcov ,STRC,thres,xyz] = myspm_graph(xSPM,SPM,hReg, cfg)
% This is a modified version of spm_graph to work with myspm_result.m
%
% [Y y beta Bcov,strc,thres,xyz] = myspm_graph(xSPM,SPM,hReg, cfg)
%
% cfg could have: 
% .plotType could be:
%           "Fitted responses" 
%           "Contrast estimates and 90% C.I."
%           "Contrast estimates and 95% C.I. effect of interest"
%
% (cc?) sgKIM, 2019.

if ~isfield(cfg,'suppressPlot'), cfg.suppressPlot=0; end

%
% Graphical display of adjusted data
%
% FORMAT [Y y beta Bcov,strc,thres,xyz] = myspm_graph(xSPM,SPM,hReg, cfg)
%
% xSPM   - structure containing SPM, distributional & filtering details
%          about the excursion set
% SPM    - structure containing generic details about the analysis
% hReg   - handle of MIP register or [x y z] coordinates
% Ic     - contrast index
% xXi    - explanatory variable index in the desginmatrix
%
%
% Y      - fitted   data for the selected voxel
% y      - adjusted data for the selected voxel
% beta   - parameter estimates (ML or MAP)
% Bcov   - Covariance of parameter estimates (ML or conditional)
% xyz    - peak coordinate in MNI-space (mm)
%
% See spm_getSPM for details.
%__________________________________________________________________________
%
% spm_graph is a Callback script that uses the structures above to:  (i)
% send adjusted (y) and fitted data (Y), for the selected voxel, to the
% workspace and (ii) provide graphics for:
%
% a) Contrasts of parameter estimates (e.g. activations) and their
% standard error.
%
% b) Fitted and adjusted responses that can be plotted against time, scan,
% or an indicator variable in the design matrix.
%
% c) (fMRI only).  Evoked responses using the basis functions to give
% impulse responses that would have been seen in the absence of other
% effects. The PSTH (peristimulus-time histogram) option provides a finite
% impulse response (FIR) estimate of the trial-specific evoked response as
% a function of peristimulus time.  This is estimated by refitting a
% convolution model to the selected voxel using an FIR basis set.  This is
% simply a set of small boxes covering successive time bins after trial
% onset.  The width of each bin is usually the TR.  This option provides a
% more time-resolved quantitative characterisation of the evoked
% hemodynamic response.  However, it should not be over-interpreted because
% inference is usually made using a simpler and more efficient basis set
% (e.g., canonical hrf, or canonical plus time derivative).
%
% Getting adjusted data:
% Ensuring the data are adjusted properly can be important (e.g. in
% constructing explanatory variables such as in a psychophysiological
% interaction). To remove or correct for specific effects, specify an
% appropriate F contrast and simply plot the fitted (and adjusted)
% responses after selecting that F contrast. The vectors Y (fitted) and y
% (adjusted) in the workspace will now be corrected for the effects in the
% reduced design matrix (X0) specified in the contrast manager with the
% column indices (iX0) of the confounds in this adjustment.
%
% Plotting data:
% All data and graphics use filtered/whitened data and residuals. In PET
% studies the parameter estimates and the fitted data are often the same
% because the explanatory variables are simply indicator variables taking
% the value of one.  Only contrasts previously defined can be plotted. This
% ensures that the parameters plotted are meaningful even when there is
% collinearity among the design matrix subpartitions.
%
% Selecting contrasts used for PPMs will automatically give plots
% based on conditonal estimates.
%
% The structure     contrast.contrast      = cbeta;
%                   contrast.standarderror = SE;
%                   contrast.interval      = 2*CI;
%
% is assigned in base workspace for plots of contrasts and their error.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% from Karl Friston's spm_graph.m
% $Id: spm_graph.m 4262 2011-03-25 13:39:20Z guillaume $

% modified by sgKIM, 2015.
Y=[]; y=[]; beta=[]; Bcov=[]; STRC=[]; thres=[]; peakxyz=[];

Ic=cfg.Ic;  % index of contrast
xXi=cfg.xXi; % column # of the variable of interest in the design matrix


%-Get Graphics figure handle
%--------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');

%-Delete previous axis and their pagination controls (if any)
%--------------------------------------------------------------------------
spm_results_ui('Clear',Fgraph,2);

%-Find nearest voxel [Euclidean distance] in point list & update GUI
%--------------------------------------------------------------------------
if isempty(xSPM.XYZmm)
  spm('alert!','No suprathreshold voxels!',mfilename,0);
  Y = []; y = []; beta = []; Bcov = [];
  return
end

% Get current coordiante (mm) from the marker
if numel(hReg) == 1
  xyz = spm_XYZreg('GetCoords',hReg);
else
  xyz = hReg;
end

% Find the nearest world coordinate from the header (xSPM.XYZmm)
[xyz,i] = spm_XYZreg('NearestXYZ', xyz, xSPM.XYZmm);
if numel(hReg) == 1, spm_XYZreg('SetCoords',xyz,hReg); end

% Find the voxel coordinate of it (xSPM.XYZ)
XYZ = xSPM.XYZ(:,i);

if isfield(cfg,'plotType')
  Cplot = cfg.plotType;
else
  Cplot = 'Fitted responses';
end

%-------------------------------
% Getting structure names
%-------------------------------
if ~isfield(cfg,'atlas')
  JOB.atlas='fsl';
end
if strcmpi(cfg.atlas,'spm12')
  [strc] = myspm_NMatlas(xyz);
  strcname = [strc.name];
  STRC=[]; STRC.strc=strc; STRC.strc_all=strc;
elseif strcmpi(cfg.atlas,'fsl')
  [strc,strc_all] = myfsl_atlasquery(xyz);
  strcname = [strc.name,' (',num2str(round(strc.prob)),'%)'];
  STRC=[]; STRC.strc=strc; STRC.strc_all=strc_all;
end


if strcmp(xSPM.STAT,'T')
  dftxt=['(',num2str(xSPM.df(end)),')'];
elseif strcmp(xSPM.STAT,'F')
  dftxt=['(',num2str(xSPM.df(1)),',',num2str(xSPM.df(end)),')'];
end
peakStat=['Peak ',xSPM.STAT,dftxt,' = ',num2str(xSPM.Z(i))];
TITLE = {strcname,peakStat};

spm('Pointer','Watch');

%-Extract filtered and whitened data from files
%==========================================================================
try
  y = spm_get_data(SPM.xY.VY, XYZ);
catch
  try
    % remap files in SPM.xY.P if SPM.xY.VY is no longer valid
    %------------------------------------------------------------------
    SPM.xY.VY = spm_vol(SPM.xY.P);
    y = spm_get_data(SPM.xY.VY, XYZ);
  end
end

if ~isempty(y), y = spm_filter(SPM.xX.K,SPM.xX.W*y); end

XYZstr = sprintf(' at [%g, %g, %g]',xyz);

%-Compute residuals
%--------------------------------------------------------------------------
if isempty(y)
  
  % make R = NaN so it will not be plotted
  %----------------------------------------------------------------------
  R   = NaN(size(SPM.xX.X,1),1);
  
else
  
  % residuals (non-whitened)
  %----------------------------------------------------------------------
  R   = spm_sp('r',SPM.xX.xKXs,y);
  
end

%-Get parameter and hyperparameter estimates
%==========================================================================
if xSPM.STAT ~= 'P'
  
  %-Parameter estimates:   beta = xX.pKX*xX.K*y;
  %-Residual mean square: ResMS = sum(R.^2)/xX.trRV
  %----------------------------------------------------------------------
  beta  = spm_get_data(SPM.Vbeta, XYZ);
  ResMS = spm_get_data(SPM.VResMS,XYZ);
  Bcov  = ResMS*SPM.xX.Bcov;
  
else
  % or conditional estimates with
  % Cov(b|y) through Taylor approximation
  %----------------------------------------------------------------------
  beta  = spm_get_data(SPM.VCbeta, XYZ);
  
  if isfield(SPM.PPM,'VB');
    % Get approximate posterior covariance at ic
    % using Taylor-series approximation
    
    % Get posterior SD beta's
    Nk=size(SPM.xX.X,2);
    for k=1:Nk,
      sd_beta(k,:) = spm_get_data(SPM.VPsd(k),XYZ);
    end
    
    % Get AR coefficients
    nsess=length(SPM.Sess);
    for ss=1:nsess,
      for p=1:SPM.PPM.AR_P
        Sess(ss).a(p,:) = spm_get_data(SPM.PPM.Sess(ss).VAR(p),XYZ);
      end
      % Get noise SD
      Sess(ss).lambda = spm_get_data(SPM.PPM.Sess(ss).VHp,XYZ);
    end
    
    % Which block are we in ?
    % this needs updating s.t xSPM contains labels of selected voxels
    v = find((SPM.xVol.XYZ(1,:)==XYZ(1))&(SPM.xVol.XYZ(2,:)==XYZ(2))&(SPM.xVol.XYZ(3,:)==XYZ(3)));
    block_index = SPM.xVol.labels(v);
    Bcov=zeros(Nk,Nk);
    for ss=1:nsess,
      % Reconstuct approximation to voxel wise correlation matrix
      post_R=SPM.PPM.Sess(ss).block(block_index).mean.R;
      if SPM.PPM.AR_P > 0
        dh=Sess(ss).a(:,1)'-SPM.PPM.Sess(ss).block(block_index).mean.a;
      else
        dh=[];
      end
      dh=[dh Sess(ss).lambda(1)-SPM.PPM.Sess(ss).block(block_index).mean.lambda];
      for i=1:length(dh),
        post_R=post_R+SPM.PPM.Sess(ss).block(block_index).mean.dR(:,:,i)*dh(i);
      end
      % Get indexes of regressors specific to this session
      scol=SPM.Sess(ss).col;
      mean_col_index=SPM.Sess(nsess).col(end)+ss;
      scol=[scol mean_col_index];
      
      % Reconstuct approximation to voxel wise covariance matrix
      Bcov(scol,scol) = Bcov(scol,scol) + (sd_beta(scol,1)*sd_beta(scol,1)').*post_R;
    end
    
  else
    Bcov     = SPM.PPM.Cby;
    for j = 1:length(SPM.PPM.l)
      l    = spm_get_data(SPM.VHp(j),XYZ);
      Bcov = Bcov + SPM.PPM.dC{j}*(l - SPM.PPM.l(j));
    end
  end
end
CI    = 1.6449;                 % = spm_invNcdf(1 - 0.05);

spm('Pointer','Arrow');
%-Plot
%==========================================================================

%-Colour specifications and index;
%--------------------------------------------------------------------------
Col   = [0 0 0; .4 .4 .4; 1 .5 .5];
switch Cplot
    %-Modeling evoked responses based on Sess
    %======================================================================
    case 'Event-related responses'
        dt    = SPM.xBF.dt;
        s     = xG.spec.Sess;
        u     = xG.spec.u;
        
        % event-related response
        %------------------------------------------------------------------
        if isempty(y)
            warning(['Data not available. ' ...
                'Plotting fitted response and 90% C.I. instead.']);
            xG.spec.Rplot = 'fitted response and 90% C.I.';
        end
        switch xG.spec.Rplot
            case 'fitted response and PSTH'
                % build a simple FIR model subpartition (X); bin size = TR
                %----------------------------------------------------------
                BIN         = SPM.xY.RT;
                xBF         = SPM.xBF;
                U           = SPM.Sess(s).U(u);
                U.u         = U.u(:,1);
                xBF.name    = 'Finite Impulse Response';
                xBF.order   = round(32/BIN);
                xBF.length  = xBF.order*BIN;
                xBF         = spm_get_bf(xBF);
                BIN         = xBF.length/xBF.order;
                X           = spm_Volterra(U,xBF.bf,1);
                k           = SPM.nscan(s);
                X           = X([0:(k - 1)]*SPM.xBF.T + SPM.xBF.T0 + 32,:);

                % place X in SPM.xX.X
                %----------------------------------------------------------
                jX          = SPM.Sess(s).row;
                iX          = SPM.Sess(s).col(SPM.Sess(s).Fc(u).i);
                iX0         = [1:size(SPM.xX.X,2)];
                iX0(iX)     = [];
                X           = [X SPM.xX.X(jX,iX0)];
                X           = SPM.xX.W(jX,jX)*X;
                X           = [X SPM.xX.K(s).X0];

                % Re-estimate to get PSTH and CI
                %----------------------------------------------------------
                j           = xBF.order;
                xX          = spm_sp('Set',X);
                pX          = spm_sp('x-',xX);
                PSTH        = pX*y(jX);
                res         = spm_sp('r',xX,y(jX));
                df          = size(X,1) - size(X,2);
                bcov        = pX*pX'*sum(res.^2)/df;
                PSTH        = PSTH(1:j)/dt;
                PST         = [1:j]*BIN - BIN/2;
                PCI         = CI*sqrt(diag(bcov(1:j,(1:j))))/dt;
        end

        % basis functions and parameters
        %------------------------------------------------------------------
        X     = SPM.xBF.bf/dt;
        x     = ([1:size(X,1)] - 1)*dt;
        j     = SPM.Sess(s).col(SPM.Sess(s).Fc(u).i(1:size(X,2)));
        B     = beta(j);

        % fitted responses with standard error
        %------------------------------------------------------------------
        Y     = X*B;
        CI    = CI*sqrt(diag(X*Bcov(j,j)*X'));

        % peristimulus times and adjusted data (y = Y + R)
        %------------------------------------------------------------------
        pst   = SPM.Sess(s).U(u).pst;
        bin   = round(pst/dt);
        q     = find((bin >= 0) & (bin < size(X,1)));
        y     = R(SPM.Sess(s).row(:));
        pst   = pst(q);
        y     = y(q) + Y(bin(q) + 1);

        % returned values
        %------------------------------------------------------------------
        if strcmp(xG.spec.Rplot,'fitted response and PSTH')
            G.PST  = PST;
            G.PSTH = PSTH;
            G.PCI  = PCI;
        end
        G.x   = x;
        G.CI  = CI;
        G.pst = pst;
  
  
%   %-Plot parameter estimates
%   %======================================================================
%   case 'Contrast estimates and 90% C.I. effect of interest'
%     Ic_eoi = find(contains({SPM.xCon.name},'Effect of interest'));
%     
%     % compute contrast of parameter estimates and 90% C.I.
%     %------------------------------------------------------------------
%     cbeta = SPM.xCon(Ic_eoi).c'*beta;
%     SE    = sqrt(diag(SPM.xCon(Ic_eoi).c'*Bcov*SPM.xCon(Ic_eoi).c));
%     CI    = CI*SE;
%     
%     contrast.contrast      = cbeta;
%     contrast.standarderror = SE;
%     contrast.interval      = 2*CI;
%     assignin('base','contrast',contrast)
%     
%     % bar chart
%     %------------------------------------------------------------------
%     figure(Fgraph)
%     subplot(2,1,2)
%     cla
%     hold on
%     
%     % estimates
%     %------------------------------------------------------------------
%     h     = bar(cbeta);
%     set(h,'FaceColor',Col(2,:))
%     
%     % standard error
%     %------------------------------------------------------------------
%     for j = 1:length(cbeta)
%       line([j j],([CI(j) -CI(j)] + cbeta(j)),...
%         'LineWidth',6,'Color',Col(3,:))
%     end
%     
%     title(TITLE,'FontSize',12)
%     xlabel('regressor')
%     ylabel(['contrast estimate',XYZstr])
%     set(gca,'XLim',[0.4 (length(cbeta) + 0.6)])
%             set(gca,'xticklabel',cat(1,SPM.Sess.U.name),...
%           'xtick',1:length(SPM.Sess.U))
%     hold off
%     
%     set(gca, 'xtick',1:numel(cbeta), 'xticklabel',cat(1,SPM.Sess(1).U.name))
%     
%     % set Y to empty so outputs are assigned
%     %------------------------------------------------------------------
%     Y = [];
  
  %-Plot parameter estimates
  %======================================================================
  case 'Contrast estimates and 90% C.I.'
    
    % compute contrast of parameter estimates and 90% C.I.
    %------------------------------------------------------------------
    cbeta = SPM.xCon(Ic).c'*beta;
    SE    = sqrt(diag(SPM.xCon(Ic).c'*Bcov*SPM.xCon(Ic).c));
    CI    = CI*SE;
    
    contrast.contrast      = cbeta;
    contrast.standarderror = SE;
    contrast.interval      = 2*CI;
    assignin('base','contrast',contrast)
    
    % bar chart
    %------------------------------------------------------------------
    figure(Fgraph)
    subplot(2,1,2)
    cla
    hold on
    
    % estimates
    %------------------------------------------------------------------
    h     = bar(cbeta);
    set(h,'FaceColor',Col(2,:))
    
    % standard error
    %------------------------------------------------------------------
    for j = 1:length(cbeta)
      line([j j],([CI(j) 0 - CI(j)] + cbeta(j)),...
        'LineWidth',6,'Color',Col(3,:))
    end
    
    title(TITLE,'FontSize',12)
    xlabel('contrast')
    ylabel(['contrast estimate',XYZstr])
    set(gca,'XLim',[0.4 (length(cbeta) + 0.6)]);
    if isfield(SPM,'Sess')
      set(gca,'xticklabel',cat(2,SPM.Sess(1).U.name),...
        'xtick',1:length(SPM.Sess(1).U))
    end
    hold off
    
    % set Y to empty so outputs are assigned
    %------------------------------------------------------------------
    Y = [];
    
    
    %-All fitted effects or selected effects
    %======================================================================
  case 'Fitted responses'
    
    % predicted or adjusted response
    %------------------------------------------------------------------
    %        str   = 'predicted or adjusted response?';
    %        if spm_input(str,'!+1','b',{'predicted','adjusted'},[1 0]);
    
    % fitted (predicted) data (Y = X1*c*c^-1*beta)
    %--------------------------------------------------------------
    Y = SPM.xX.X*SPM.xCon(Ic).c*pinv(SPM.xCon(Ic).c)*beta;
    
    % adjusted data
    %------------------------------------------------------------------
    y     = Y + R;
        
    % when contrast has one only non-zero value:
    i    = xXi;
    %x    = SPM.xX.xKXs.X(:,xXi);
    x  = SPM.xX.xKXs.X * cfg.contrast;
    XLAB = SPM.xX.name{i};
    if ~exist('skipTheGraph','var')
      %offset adjustment (by sgKIM)
      yoffset = beta(1);
      y0 = y + yoffset;
      Y1 = Y + yoffset;
      try xoffset = mean(evalin('base','JOB.vi.val'));
      catch ME
        xoffset=0;
      end
      
      x = x + xoffset;
      
      if ~cfg.suppressPlot
        
        % plot
        %------------------------------------------------------------------
        figure(Fgraph)
        subplot(2,1,2)
        cla
        hold on
        [p q] = sort(x);
        %Col(4,:)=[.5 .5 1];
        h=[0 0];
        if mod(cfg.Ic,2) == 0
          x=-x;
        end
        
        if cfg.isfmri
          xG=[];
          xG.def='Fitted responses';
          xG.spec.Ic = Ic;
          xG.spec.predicted = 1;
          xG.spec.x.scan = 1;
          [Y,y,beta,Bcov,G] = spm_graph(SPM,XYZ,xG);
          h(1)=plot(G.x,zscore(y),'.-','MarkerSize',8, 'Color',Col(3,:)); % observed
          h(2)=plot(G.x,zscore(Y),'LineWidth',2,'Color',Col(2,:));        % predicted
          cfg.x_name='Time [s]';
          cfg.y_name='Z(y)';
          
        elseif all(diff(x(q))) % no duplication of x (thus likely continuous..?)
          h(1)=plot(x(q),y0(q),'o','MarkerSize',8, 'Color',Col(3,:)); %offset adjustment (by sgKIM)
          h(2)=plot(x(q),Y1(q),'LineWidth',2,'Color',Col(2,:));
        
        else % for discrete values
          try h(2)=plot(x(q),Y1(q),'o','MarkerSize',8,'Color',Col(1,:));
          catch ME
            hh=get(gca,'children');
            h(2)=hh(end);
          end
          %h(2)=plot(x(q),y(q),'.','MarkerSize',8, 'Color',Col(2,:));
          x_fit=x(q);
          y_fit=y0(q);
          [~,idx]=sort(x_fit);
          try h(1)=plot(x(q),y0(q),'o','MarkerSize',10, 'Color',Col(3,:)); %offset adjustment (by sgKIM)
          catch ME
            hh=get(gca,'children');
            h(1)=hh(end);
          end
          xlim = get(gca,'XLim');
          xlim = [-1 1]*diff(xlim)/4 + xlim;
          set(gca,'XLim',xlim)
        end
        
        plotdata=[];
        plotdata.x = x(q);
        plotdata.y = y0(q);
        plotdata.yhat = Y1(q);
        if cfg.isfmri
          plotdata.x = G.x;
          plotdata.y = y;
          plotdata.yhat = Y;
        end
        save ([xSPM.swd,'/plotdata_cluster',num2str(cfg.clusterIdx),'.mat'],'plotdata');
        
        
        if isfield(cfg,'markCorrThres')
          % ref: http://en.wikipedia.org/wiki/Bonferroni_correction
          % Z=atanh(R) ~ N(0,1/sqrt(n-3)) where n is # of samples
          % As we compute the correlation from 410 volumes...
          
          % This is bonferoni correction:
          if isfield(cfg.markCorrThres, 'CorrFDRThres')
            thres=cfg.markCorrThres.CorrFDRThres;
          else
            % How about FDR correction?
            alphalevel=cfg.markCorrThres.alpha;
            NumFrames=cfg.markCorrThres.NumFrames;
            SE=1/sqrt(NumFrames-3);
            % that means, we set a threshold for each individual
            NumSubj=size(SPM.xY.P,1);
            thres=zeros(NumSubj,1);
            for subji=1:NumSubj
              r=sort(spm_get_data(SPM.xY.VY(subji), SPM.xVol.XYZ));
              z=atanh(r);
              pvals = min(normcdf(z,0,SE), normcdf(-z,0,SE))*2;
              thres(subji) = r(find(pvals > fdr(pvals,alphalevel),1,'last'));
            end
            thres=mean(thres);
          end
          CA=axis;
          h_alpha=line([CA(1:2); CA(1:2)]', [thres thres; -thres -thres]', 'color',[.2 .8 .2]);
        else
          thres=[];
        end
        
        title(TITLE,'FontSize',12)
        
        if isfield(cfg,'x_name')
          xlabel(cfg.x_name);
        else
          xlabel(XLAB)
        end
        if isfield(cfg,'y_name')
          ylabel([cfg.y_name,XYZstr])
        else
          ylabel(['response',XYZstr])
        end
        [~,a]=min(x);
        [~,b]=max(x);
        if y(a) < y(b)
          LOC='southeast';
        else
          LOC='northeast';
        end
        if exist('h_alpha','var')
          legend([h h_alpha'], {'observed','fitted',['q=',sprintf('%0.2f',cfg.markCorrThres.alpha)]},'location',LOC);
        else
          legend(h, {'observed','fitted'},'location',LOC);
        end
        hold off
        if cfg.isfmri
          % for fMRI, show only first <200 TRs (otherwise it's too dense to see)
          % FIND the first peak of the regressor
          idx1 = max([1 find(Y~=0,1,'first')-20]);
          idx2 = min([numel(plotdata.x) 180+idx1-1]);
          % SET xlim for first peak, and 100 TRs later:
          set(gca,'xlim',[plotdata.x(idx1) plotdata.x(idx2)])
          grid on; box on;
        end
      end

      %%
    end
end

% Turn hold button off - this will alert the user to press it again
%--------------------------------------------------------------------------
try
  set(get(gcbo,'Userdata'),'Value',0);
catch
end


%-call Plot UI
%--------------------------------------------------------------------------
spm_results_ui('PlotUi',gca)
end

%% myfsl_atlasquery.m

function [strc, strc_all]=myfsl_atlasquery(mni_xyz, ijkflag, atlasset)
% [strc, strc_all]=myfsl_atlasquery(mni_xyz, ijkflag)
% mni_xyz: MNI-coordinates (mm) or 1-based ijk (voxels) with nifti
%
% will find name for a given MNI coordiante (x,y,z) mm
%   or voxel index (i,j,k) by marking ijkflag
% with no input arguments, will try to read a coordinate from SPM figure.
%
% (cc) 2015, 2019. sgKIM. solleo@gmail.com

if ~exist('ijkflag','var')
  ijkflag=0;
end

if ~exist('atlasset','var')
  atlasset='gm';
end

fslpath = getenv('FSLDIR');
xmlFNAMES{1}=fullfile(fslpath,'data/atlases/HarvardOxford-Cortical-Lateralized.xml');
niiFNAMES{1}=fullfile(fslpath,'data/atlases/HarvardOxford/HarvardOxford-cortl-prob-2mm.nii.gz');
if ~exist(niiFNAMES{1},'file')
  xmlFNAMES{1}=fullfile(fslpath,'data/atlases/HarvardOxford-Cortical.xml');
  niiFNAMES{1}=fullfile(fslpath,'data/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz');
end
xmlFNAMES{2}=fullfile(fslpath,'data/atlases/HarvardOxford-Subcortical.xml');
niiFNAMES{2}=fullfile(fslpath,'data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz');
xmlFNAMES{3}=fullfile(fslpath,'data/atlases/Cerebellum_MNIfnirt.xml');
niiFNAMES{3}=fullfile(fslpath,'data/atlases/Cerebellum/Cerebellum-MNIfnirt-prob-2mm.nii.gz');
xmlFNAMES{4}=fullfile(fslpath,'data/atlases/JHU-tracts.xml');
niiFNAMES{4}=fullfile(fslpath,'data/atlases/JHU/JHU-ICBM-tracts-prob-2mm.nii.gz');

% okay unzipping takes more than 2 sec..., so try to find it from /tmp first
for a=1:4
  [~,fname1,~]=fileparts(niiFNAMES{a});
  if ~exist(['/tmp/',fname1], 'file')
    gunzip(niiFNAMES{a}, '/tmp/');
  end
  niiFNAMES{a} = ['/tmp/',fname1];
  xDOC{a}=xml2struct(xmlFNAMES{a});
end

%% get ijk coordinate
if ~ijkflag
  if nargin == 0
    hReg= evalin('base','hReg;');
    xSPM= evalin('base','xSPM;');
    if numel(hReg) == 1
      xyz = spm_XYZreg('GetCoords',hReg);
    else
      xyz = hReg;
    end
  else
    xyz = mni_xyz';
  end
  try xSPM.XYZmm
    [xyz,i] = spm_XYZreg('NearestXYZ', xyz ,xSPM.XYZmm);
  catch ME
    xyz = mni_xyz';
  end
  ijk = round(xyz2ijk(xyz, niiFNAMES{1}))';
else
  ijk = mni_xyz';
  xyz = round(ijk2xyz(ijk', niiFNAMES{1}))';
end

%% now read probs for XYZ using spm_get_data (very efficient when reading only one voxel)
strc.name='N/A';
strc.prob=0;
strc_all.name=cell(1,4);
strc_all.prob=[0 0 0 0];
k=1;
if strcmpi('all',atlasset)
  ATLAS=1:4;
elseif strcmpi('gm',atlasset)
  ATLAS=1:3;
elseif strcmpi('wm',atlasset)
  ATLAS=4;
end
for a=ATLAS % for each atlas
  P = spm_vol(niiFNAMES{a});
  probs = spm_get_data(P, ijk);
  if a==2
    % ignore cerebral grey/white matter in the subcortical atlas
    probs([1 2 12 13],:)=0;
  end
  
  % find maximal prob from
  if size(ijk,2) == 1 % this works for probs<Nx1>
    [~,b]=max(probs);
    probs_b = probs(b);
  else % in case of cluster.. probs<NxV>
    % I want to compute mean prob for each label
    probs_ri = round(mean(probs,2));
    [~,b]=max(probs_ri);
    probs_b = probs_ri(b);
  end
  
  if probs_b > strc.prob % is the current maximal greater than the previous one?
    strc.name = xDOC{a}.atlas.data.label{b}.Text;
    strc.prob = probs_b;
  end
  
  % for any other possibilities...
  if size(ijk,2) == 1 % for a voxel
    nz = find(~~probs); % non-zero probs.
    for j=1:numel(nz)
      strc_all.name{k}=xDOC{a}.atlas.data.label{nz(j)}.Text;
      strc_all.prob(k)=probs(nz(j));
      k=k+1;
    end
  else % for a cluster
    nz = find(~~probs_ri);
    for j=1:numel(nz)
      strc_all.name{k}=xDOC{a}.atlas.data.label{nz(j)}.Text;
      strc_all.prob(k)=probs_ri(nz(j));
      k=k+1;
    end
  end
end

strc_unsort=strc_all;
% and sort;
[~,idx] = sort(strc_unsort.prob, 'descend');
for j=1:numel(idx)
  strc_all.name{j} = strc_unsort.name{idx(j)};
  strc_all.prob(j) = strc_unsort.prob(idx(j));
end

end


%% =============================================================================
% source: http://www.mathworks.com/matlabcentral/fileexchange/28518-xml2struct

function [ s ] = xml2struct( file )
%Convert xml file into a MATLAB structure
% [ s ] = xml2struct( file )
%
% A file containing:
% <XMLname attrib1="Some value">
%   <Element>Some text</Element>
%   <DifferentElement attrib2="2">Some more text</Element>
%   <DifferentElement attrib3="2" attrib4="1">Even more text</DifferentElement>
% </XMLname>
%
% Will produce:
% s.XMLname.Attributes.attrib1 = "Some value";
% s.XMLname.Element.Text = "Some text";
% s.XMLname.DifferentElement{1}.Attributes.attrib2 = "2";
% s.XMLname.DifferentElement{1}.Text = "Some more text";
% s.XMLname.DifferentElement{2}.Attributes.attrib3 = "2";
% s.XMLname.DifferentElement{2}.Attributes.attrib4 = "1";
% s.XMLname.DifferentElement{2}.Text = "Even more text";
%
% Please note that the following characters are substituted
% '-' by '_dash_', ':' by '_colon_' and '.' by '_dot_'
%
% Written by W. Falkena, ASTI, TUDelft, 21-08-2010
% Attribute parsing speed increased by 40% by A. Wanner, 14-6-2011
% Added CDATA support by I. Smirnov, 20-3-2012
%
% Modified by X. Mo, University of Wisconsin, 12-5-2012

if (nargin < 1)
  clc;
  help xml2struct
  return
end

if isa(file, 'org.apache.xerces.dom.DeferredDocumentImpl') || isa(file, 'org.apache.xerces.dom.DeferredElementImpl')
  % input is a java xml object
  xDoc = file;
else
  %check for existance
  if (exist(file,'file') == 0)
    %Perhaps the xml extension was omitted from the file name. Add the
    %extension and try again.
    if (isempty(strfind(file,'.xml')))
      file = [file '.xml'];
    end
    
    if (exist(file,'file') == 0)
      error(['The file ' file ' could not be found']);
    end
  end
  %read the xml file
  xDoc = xmlread(file);
end

%parse xDoc into a MATLAB structure
s = parseChildNodes(xDoc);

end

% ----- Subfunction parseChildNodes -----
function [children,ptext,textflag] = parseChildNodes(theNode)
% Recurse over node children.
children = struct;
ptext = struct; textflag = 'Text';
if hasChildNodes(theNode)
  childNodes = getChildNodes(theNode);
  numChildNodes = getLength(childNodes);
  
  for count = 1:numChildNodes
    theChild = item(childNodes,count-1);
    [text,name,attr,childs,textflag] = getNodeData(theChild);
    
    if (~strcmp(name,'#text') && ~strcmp(name,'#comment') && ~strcmp(name,'#cdata_dash_section'))
      %XML allows the same elements to be defined multiple times,
      %put each in a different cell
      if (isfield(children,name))
        if (~iscell(children.(name)))
          %put existsing element into cell format
          children.(name) = {children.(name)};
        end
        index = length(children.(name))+1;
        %add new element
        children.(name){index} = childs;
        if(~isempty(fieldnames(text)))
          children.(name){index} = text;
        end
        if(~isempty(attr))
          children.(name){index}.('Attributes') = attr;
        end
      else
        %add previously unknown (new) element to the structure
        children.(name) = childs;
        if(~isempty(text) && ~isempty(fieldnames(text)))
          children.(name) = text;
        end
        if(~isempty(attr))
          children.(name).('Attributes') = attr;
        end
      end
    else
      ptextflag = 'Text';
      if (strcmp(name, '#cdata_dash_section'))
        ptextflag = 'CDATA';
      elseif (strcmp(name, '#comment'))
        ptextflag = 'Comment';
      end
      
      %this is the text in an element (i.e., the parentNode)
      if (~isempty(regexprep(text.(textflag),'[\s]*','')))
        if (~isfield(ptext,ptextflag) || isempty(ptext.(ptextflag)))
          ptext.(ptextflag) = text.(textflag);
        else
          %what to do when element data is as follows:
          %<element>Text <!--Comment--> More text</element>
          
          %put the text in different cells:
          % if (~iscell(ptext)) ptext = {ptext}; end
          % ptext{length(ptext)+1} = text;
          
          %just append the text
          ptext.(ptextflag) = [ptext.(ptextflag) text.(textflag)];
        end
      end
    end
    
  end
end
end

% ----- Subfunction getNodeData -----
function [text,name,attr,childs,textflag] = getNodeData(theNode)
% Create structure of node info.

%make sure name is allowed as structure name
name = toCharArray(getNodeName(theNode))';
name = strrep(name, '-', '_dash_');
name = strrep(name, ':', '_colon_');
name = strrep(name, '.', '_dot_');

attr = parseAttributes(theNode);
if (isempty(fieldnames(attr)))
  attr = [];
end

%parse child nodes
[childs,text,textflag] = parseChildNodes(theNode);

if (isempty(fieldnames(childs)) && isempty(fieldnames(text)))
  %get the data of any childless nodes
  % faster than if any(strcmp(methods(theNode), 'getData'))
  % no need to try-catch (?)
  % faster than text = char(getData(theNode));
  text.(textflag) = toCharArray(getTextContent(theNode))';
end

end

% ----- Subfunction parseAttributes -----
function attributes = parseAttributes(theNode)
% Create attributes structure.

attributes = struct;
if hasAttributes(theNode)
  theAttributes = getAttributes(theNode);
  numAttributes = getLength(theAttributes);
  
  for count = 1:numAttributes
    %attrib = item(theAttributes,count-1);
    %attr_name = regexprep(char(getName(attrib)),'[-:.]','_');
    %attributes.(attr_name) = char(getValue(attrib));
    
    %Suggestion of Adrian Wanner
    str = toCharArray(toString(item(theAttributes,count-1)))';
    k = strfind(str,'=');
    attr_name = str(1:(k(1)-1));
    attr_name = strrep(attr_name, '-', '_dash_');
    attr_name = strrep(attr_name, ':', '_colon_');
    attr_name = strrep(attr_name, '.', '_dot_');
    attributes.(attr_name) = str((k(1)+2):(end-1));
  end
end
end

%% ------ SG's functions
% 
% function xyz = ijk2xyz(ijk, nii)
% % xyz = ijk2xyz(ijk, nii)
% % converts world-to-voxel coordinates.
% %
% % Inputs:
% %   ijk   [Nx3 vector] is 1-based voxel coordinates for MATLAB
% % Output:
% %   xyz   [Nx3 vector] is world-coordinates (e.g. MNI-coord)
% %   nii   the nii structure read using load_untouch_nii.m (or the filename)
% %
% % see xyz2ijk.m
% % (cc) sgKIM, 2014. solleo@gmail.com
% 
% 
% if ischar(nii)
%   nii = load_uns_nii(nii);
% end
% 
% T=[nii.hdr.hist.srow_x; nii.hdr.hist.srow_y; nii.hdr.hist.srow_z; 0 0 0 1];
% 
% if numel(ijk) == 3
%   if size(ijk,1) > size(ijk,2)
%     ijk=ijk';
%   end
% else
%   if size(ijk,1) < size(ijk,2)
%     ijk=ijk';
%   end
% end
% 
% xyz = (T)*[(ijk-1) ones(size(ijk,1),1) ]';
% xyz(4,:)=[];
% xyz=xyz';
% 
% end
% 
% %%  ijk = xyz2ijk(xyz, nii)
% % converts world-to-voxel coordinates.
% %
% % Inputs:
% %   xyz   [Nx3 vector] is world-coordinates (e.g. MNI-coord)
% %   nii   the nii structure read using load_untouch_nii.m (or the filename)
% % Output:
% %   ijk   [Nx3 vector] is 1-based voxel coordinates for MATLAB
% %
% % see ijk2xyz.m
% % (cc) sgKIM, 2014. solleo@gmail.com
% 
% % function ijk = xyz2ijk(xyz, nii)
% % 
% % if ischar(nii)
% %   hdr= load_untouch_header_only(nii);
% % else
% %   hdr = nii.hdr;
% % end
% % 
% % if numel(xyz) == 3
% %   if size(xyz,1) > size(xyz,2)
% %     xyz=xyz';
% %   end
% % else
% %   if size(xyz,1) < size(xyz,2)
% %     xyz=xyz';
% %   end
% % end
% % 
% % T=[hdr.hist.srow_x; hdr.hist.srow_y; hdr.hist.srow_z; 0 0 0 1];
% % ijk = inv(T) * [xyz ones(size(xyz,1),1)]';
% % ijk(4,:)=[];
% % ijk = ijk'+1;
% % 
% % end
% 
% 
% function Ijk1 = xyz2ijk(Xyz, nii)
% % converts world-coordinates (mm) to 1-based voxel-subscripts
% % 
% % SYNTAX
% % Ijk1 = xyz2ijk(Xyz, Info)
% %
% % INPUTS
% % Xyz    [3xN]  World-coordinates (mm) of N points
% % Info   (1x1)  a structure from NIFTIINFO
% %      | '1xN'  a filename of the NIFTI file
% %
% % OUTPUT
% % Ijk1   [3xN]  MATLAB's 1-based voxel-subscripts of N points
% %
% % (cc) 2022, dr.seunggoo.kim@gmail.com
% 
% Info = niftiinfo(nii);
% if ischar(Info)
%   Info = niftiinfo(Info);
% else
%   error('Info can be either a structure or a string')
% end
% Ijk0 = [Xyz ones(size(Xyz,1),1)]/Info.Transform.T;
% Ijk1 = Ijk0(:,1:3) + 1;
% 
% end


%%
function TEST()
FnMni = fullfile(getenv('FSLDIR'),'data','standard','MNI152_T1_1mm.nii.gz');
ijk1 = xyz2ijk([0 0 0], FnMni); % world origin
assert(isequal(ijk1, [90 126 72]+1))
ijk1 = xyz2ijk([90 -126 -72], FnMni)
assert(isequal(ijk1, [0 0 0]+1)); % voxel origin
end

