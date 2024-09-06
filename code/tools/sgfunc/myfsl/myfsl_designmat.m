function JOB = myfsl_designmat(JOB)
% JOB = myfsl_designmat(JOB)
%
% creates .con and .mat files for randomization using FSL
%
% JOB requires:
%  .model
%  .cidx
% (.contrast)
% (.contrastname)
% (.model_desc)
%
% (cc) 2015, sgKIM.   solleo@gmail.com   https://ggooo.wordpress.com


% find regressor names
varnames = char(JOB.model);
X = double(JOB.model);
[numSubj,numReg] = size(X);

% where is the regressor to test?
if isfield(JOB,'cidx')
if numel(JOB.cidx) == 1
cidx = sort(JOB.cidx);
JOB.vi=varnames{cidx};
elseif numel(JOB.cidx) == 2
cidx = JOB.cidx;
JOB.vi = [varnames{cidx(1)},'-',varnames{cidx(2)}];
end
else
cidx = numReg;
warning(['# contrast not given explicitly, assuming ',varnames{cidx},' to test']);
JOB.vi=varnames{cidx};
end

% demeaning regressors (except binary variables and constant)
for j=1:numReg
if numel(unique(X(:,j)))>2        % if not constant or binary
X(:,j) = X(:,j) - mean(X(:,j)); % demean it
end
end

% PP height? not for TBSS but for FEAT
PPheights = range(X);

% create contrast (if not given)
if ~isfield(JOB,'contrast')
if numel(cidx) == 1
K=[
zeros(1,cidx-1),  1, zeros(1,numReg-cidx)
zeros(1,cidx-1), -1, zeros(1,numReg-cidx)
];
elseif numel(cidx) == 2
K=[
zeros(1,cidx(1)-1), +1, zeros(1,cidx(2)-cidx(1)-1), -1, zeros(1,numReg-cidx(2))
zeros(1,cidx(1)-1), -1, zeros(1,cidx(2)-cidx(1)-1), +1, zeros(1,numReg-cidx(2))
];
end
else
K = JOB.contrast;
end

% create contrast name (if not given)
if ~isfield(JOB,'contrastname')
if numel(cidx) == 1
JOB.contrastname =  {['+',char(varnames{cidx})],['-',char(varnames{cidx})]};
elseif numel(cidx) == 2
JOB.contrastname = {
[varnames{cidx(1)},'-gt-',varnames{cidx(2)}], ...
[varnames{cidx(2)},'-gt-',varnames{cidx(1)}]};
end
end

%% write .mat
if ~isfield(JOB,'model_desc'),
JOB.model_desc = fsss_model_desc(JOB.model, JOB.cidx);
end

filename = [JOB.dir_glm,'/',JOB.model_desc,'.mat'];
JOB.fname_mat = filename;
fid = fopen(filename,'w');
fprintf(fid,'/NumWaves\t%d\n',  numReg);
fprintf(fid,'/NumPoints\t%d\n', numSubj);
fprintf(fid,'/PPheights\t');
fclose(fid);
dlmwrite(filename, PPheights, 'delimiter','\t','-append');
fid = fopen(filename,'a');
fprintf(fid,'/Matrix\n');
fclose(fid);
dlmwrite(filename, X, 'delimiter','\t','-append','precision','%.8f');

%% wirte .con
filename= [JOB.dir_glm,'/',JOB.model_desc,'.con'];
JOB.fname_con = filename;
fid = fopen(filename,'w');
for j=1:numel(JOB.contrastname)
fprintf(fid,'/ContrastName%d\t"%s"\n', j, JOB.contrastname{j});
end
fprintf(fid,'/NumWaves\t%d\n', numReg);
fprintf(fid,'/NumContrasts\t%d\n', numel(JOB.contrastname));
fprintf(fid,'/Matrix\n');
fclose(fid);
dlmwrite(filename, K, 'delimiter','\t','-append','precision','%.8f');
end
