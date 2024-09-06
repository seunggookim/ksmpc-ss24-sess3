function JOB = myspm_cntrst (JOB)
% creates CON and SPM images, to be used by MYSPM_GLM and MYSPM_FMRIGLM
% and calls MYSPM_RESULT
%
% JOB = myspm_cntrst (JOB)
%
%
% JOB requires:
%  .dir_glm
%  .cntrstMtx  (for T-contrasts)
%  .titlestr   (for T-contrasts)
%  .effectOfInterest [1x1] 0 | 1 (default=1)
%  .FcntrstMtx (for F-contrasts)
%  .Ftitlestr  (for F-contrasts)
%
% (cc) 2018, 2019, sgKIM. solleo@gmail.com


spm('defaults','fmri');
if ~isfield(JOB,'effectOfInterest')
  JOB.effectOfInterest = 1;
end
load([JOB.dir_glm,'/SPM.mat'],'SPM')
if ~isfield(JOB,'isfmri')
  JOB.isfmri = isfield(SPM,'Sess');
end

%% SET DEFAULT CONTRASTS (all regressors):
if ~isfield(JOB,'cntrstMtx')
  if JOB.isfmri % 1-level fmri GLM:
    nRegInt = [];
    for j = 1:numel(SPM.Sess)
      nRegInt(j) = numel(SPM.Sess(j).U);
    end
    if any(nRegInt~=nRegInt(1))
      error(['# of regressors of interest are various across sessions!'])
    end
    k = length(SPM.Sess(1).U);
    JOB.cntrstMtx = [ones(1,k); -ones(1,k);
      kron(eye(k),[1 -1]')];
    if ~isfield(JOB,'titlestr')
      varnames = [SPM.Sess(1).U(:).name];
      JOB.titlestr = {'+All','-All'};
      for i=1:numel(varnames)
        JOB.titlestr = [JOB.titlestr, ['+',varnames{i}]];
        JOB.titlestr = [JOB.titlestr, ['-',varnames{i}]];
      end
    end
  else % 2-level GLM:
    [n, k] = size(SPM.xX.X);
    JOB.cntrstMtx = kron([zeros(k-1,1) eye(k-1)],[1 -1]');
    if ~isfield(JOB,'titlestr')
      if ~isempty(SPM.xC)
        varnames = {SPM.xC.rcname};
      else
        varnames= {'1'};
      end
      JOB.titlestr = {};
      for i=1:numel(varnames)
        JOB.titlestr = [JOB.titlestr, ['+',varnames{i}]];
        JOB.titlestr = [JOB.titlestr, ['-',varnames{i}]];
      end
    end
  end
end

%% SET T-contrasts:
matlabbatch={};
con=[];
con.spmmat = {fullfile(JOB.dir_glm, 'SPM.mat')};
if ~isfield(JOB,'NumSess')
  try
    JOB.NumSess = numel(SPM.Sess);
  catch
    JOB.NumSess = 1;
  end
end
NumSess = JOB.NumSess;
if isfield(JOB,'cntrstMtx')
  NumCnt = size(JOB.cntrstMtx,1);
  for k=1:NumCnt
    if isfield(JOB,'titlestr')
      con.consess{k}.tcon.name = JOB.titlestr{k};
    else
      con.consess{k}.tcon.name = ['Contrast#',num2str(k)];
    end
    con.consess{k}.tcon.convec = JOB.cntrstMtx(k,:);
    if NumSess>1
      con.consess{k}.tcon.sessrep = 'repl';
    else
      con.consess{k}.tcon.sessrep = 'none';
    end
    if isfield(JOB,'sessrep')
      con.consess{k}.tcon.sessrep=JOB.sessrep;
    end
  end
end


%% Effect of interest
if JOB.effectOfInterest
  if ~isfield(JOB,'FcntrstMtx')
    JOB.FcntrstMtx={};
    JOB.Ftitlestr={};
  end
  JOB.Ftitlestr = [JOB.Ftitlestr 'Effect of interest'];
  if JOB.isfmri
    for j = 1:numel(SPM.Sess)
      nRegInt(j) = numel(SPM.Sess(j).U); % number of regressors of interest
      nxBForder = size(SPM.xBF.bf,2); % order of the response function
%       nRegNsn(j) = size(SPM.Sess(j).C.C,2); % number of regressors of nuisance
    end
%     % SANITY CHECK:
%     if sum(nRegInt)*nxBForder+sum(nRegNsn)+numel(SPM.Sess) ~= size(SPM.xX.X,2)
%       error('# of regressors does not fit! CANNOT set up Effect of Interest contrast!')
%     end
%     eoicont = zeros(nRegInt(j)*nxBForder, size(SPM.xX.X,2));
%     idxRows = [0 cumsum(nRegInt*nxBForder)];
%     idxCols = [0 cumsum(nRegInt*nxBForder + nRegNsn)];
%     for j = 1:numel(SPM.Sess)
%       eoicont(1+idxRows(j):idxRows(j+1), 1+idxCols(j):idxCols(j+1)) ...
%         = [eye(nRegInt(j)*nxBForder) zeros(nRegInt(j)*nxBForder,nRegNsn(j))];
%     end
    eoicont = eye(nRegInt(j)*nxBForder);
    JOB.FcntrstMtx = [JOB.FcntrstMtx; eoicont];
  else % 2nd-level or higher
    warning('THINK ABOUT EFFECT of INTEREST for 2-level ANALYSIS is USEFUL!')
  end
end


%% SET F-contrasts
if isfield(JOB,'FcntrstMtx')
  k = NumCnt;
  for j=1:numel(JOB.FcntrstMtx)
    k = k + 1;
    con.consess{k}.fcon.name = JOB.Ftitlestr{j};
    con.consess{k}.fcon.convec = JOB.FcntrstMtx{j};
    if NumSess>1
      con.consess{k}.fcon.sessrep = 'repl';
    else
      con.consess{k}.fcon.sessrep = 'none';
    end
    if isfield(JOB,'sessrep')
      con.consess{k}.fcon.sessrep=JOB.sessrep;
    end
  end
end


%% CLEAR previous results:
if isfield(JOB,'newContrast')
  con.delete = JOB.newContrast;
else
  con.delete = 1;
end
if con.delete == 1
  unix(['rm -f ',JOB.dir_glm,'/con*']);
  unix(['rm -f ',JOB.dir_glm,'/sigclus*']);
end
matlabbatch{1}.spm.stats.con = con;


%% RUN:
spm_jobman('initcfg')
spm_jobman('run', matlabbatch)


%% Now create result reports
JOB.mygraph.y_name = 'y';
JOB.mygraph.x_name = 'x';
if ~isfield(JOB,'thres')
  JOB.thres.desc  = 'cluster';
  JOB.thres.alpha = 0.05;
end
if ~isfield(JOB,'NOREPORT') && ~isfield(JOB,'noreport')
  JOB = myspm_result(JOB);
else
%   myps2pdf(JOB.fname_spm_fig)
end


end
