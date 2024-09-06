function JOB = myspm_strcNterm(JOB)
% converts 'model' structure from SurfStat into the SPM's structure
%
% (cc) 2016. sgKIM.  mailto://solleo@gmail.com  https://ggooo.wordpress.com/

if isfield(JOB,'model')&&~isfield(JOB,'vi')
  % terms to spm.struct
  varnames = char(JOB.model);
  X = double(JOB.model);
  [numSubj,numReg] = size(X);
  % where is the regressor to test?
  if isfield(JOB,'cidx')
    c = JOB.cidx;
  else
    c = numReg;
  end
  if isempty(varnames{c})
    JOB.vi.name=['var',num2str(c)];
  else
    JOB.vi.name = varnames{c};
  end
  JOB.vi.val  = X(:,c);
  % where is the constant?
  d = find(sum(X)==numSubj);
  % what's left for covariates?
  idx = 1:numReg;
  idx([c,d]) = [];
  if numel(idx)
    for i=1:numel(idx)
      j=idx(i);
      JOB.vn(i).name = varnames{j};
      JOB.vn(i).val  = X(:,j);
    end
  end
elseif ~isfield(JOB,'model')&&isfield(JOB,'vi')
  % spm.struct to terms
  JOB.model = 1 + term(JOB.vi.val, JOB.vi.name);
  JOB.cidx = 2;
  if isfield(JOB,'vn')
    for j=1:numel(JOB.vn)
      JOB.model = JOB.model + term(JOB.vn(j).val, JOB.vn(j).name);
    end
  end
end
if ~isfield(JOB,'model_desc')
  if isfield(JOB,'cidx')
    JOB.model_desc = fsss_model_desc(JOB.model, JOB.cidx);
  else
    JOB.model_desc = fsss_model_desc(JOB.model);
  end
end
end


function [model_desc,contrastOfInterest]  = fsss_model_desc (model, cidx)
% model_desc  = fsss_model_desc (model)
%
% generates string describing a given model in "term" format
% (cc) 2015, sgKIM.   solleo@gmail.com

if isnumeric(model)
  if model == 1
    model_desc='1';
  else
    error('Use terms except the interceptor')
  end
else
  N = numel(model);
  if N == 1  % T-stat with model
    d = size(model);
    regname = char(model);
    model_desc=[];
    if numel(regname)>0
      for i=1:numel(regname)
        if isempty(regname{i})
          regname{i} = ['x',num2str(i)];
        end
        model_desc = [model_desc,'+',regname{i}];
      end
      model_desc(1)=[];
    end
    if exist('cidx','var')
      if numel(cidx) == 1
        contrastOfInterest = regname{cidx};
      elseif numel(cidx) == 2
        contrastOfInterest = [regname{cidx(1)},'-',regname{cidx(2)}];
      end
      model_desc = [model_desc,'_',contrastOfInterest];
    end
  elseif N  == 2  % F-stat with model{1} and model{2}
    d1 = size(model{1});
    d2 = size(model{2});
    if d1(2) < d2(2)
      M = model{2};
      model{2} = model{1};
      model{1} = M;
      clear M;
    end
    str1 = char(model{1});
    str2 = char(model{2});
    model_desc=[];
    for i=1:numel(str2)
      if isempty(str2{i})
        str2{i} = ['x',num2str(i-1)];
      end
      model_desc = [model_desc,'+',str2{i}];
    end
    model_desc = [model_desc,'+[',str1{i+1},']'];
  end
end

end