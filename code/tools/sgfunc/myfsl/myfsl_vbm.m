function JOB = myfsl_vbm (JOB)
% JOB requires:
% .files_query
% .meastype
% .template


%% 0. set defaults and parse input
setenv('FSLOUTPUTTYPE','NIFTI_GZ');

% 0-1. Set directories
if isfield(JOB,'dir_output')
dir_output=JOB.dir_output;
else
dir_output=pwd;
end

dir_log = [dir_output,'/job_logs'];
[~,~]=mkdir(dir_log);
cd(dir_log);

dir_figs  = [dir_output ,'/figs'];
[~,~]=mkdir(dir_figs);


%% 1-1. find input files
if ~isfield(JOB,'meastype')
error('You need to enter measurement type: "FA", "MD", "L1", "L2", "L3", or "GMD"')
else
meastype = JOB.meastype;
end

if isfield(JOB,'files_query')
[~,fnames] = mydir(JOB.files_query);
elseif isfield(JOB,'filenames')
fnames = JOB.filenames;
else
error('You need to specify inputs in JOB.files_query or JOB.filenames');
end
disp('> Input filenames:');
for j=1:numel(fnames)
ls(fnames{j});
end

% 1-2. find template
if ~isfield(JOB,'fname_template')
disp('# No template defined:');
if ~strcmpi(meastype,'gmd')
fname_template=[getenv('FSLDIR'),'/data/standard/FMRIB58_FA_1mm.nii.gz'];
else
fname_template=[getenv('FSLDIR'),'/data/standard/MNI152_T1_1mm_brain.nii.gz'];
end
else
fname_template = JOB.fname_template;
end
disp('> Set template:');
ls(fname_template);


%% 2. create design matrix and contrasts

% 2-1. find cidx: (index of contrast variable)
if ~isfield(JOB,'contrasts')
if isfield(JOB,'cidx')
cidx = JOB.cidx;
end
else
cidx = find(JOB.contrasts(1,:));
end

% 2-2. make directory
if isfield(JOB,'fwhm_mm')
fwhmstr = ['_fwhm,',num2str(JOB.fwhm_mm),'mm'];
else
fwhmstr =[];
end

dir_glm = [dir_output,'/GLM/',meastype,'_',fsss_model_desc(JOB.model, cidx),fwhmstr,'/'];
[~,~] = mkdir(dir_glm);
dir0=pwd;
cd(dir_glm);

% 2-3. create design matrix and contrast vectors
NR=size(JOB.model,2); % # of regressors
varnames = char(JOB.model);
JOB.contrasts = [zeros(1,cidx-1),  1, zeros(1,NR-cidx); zeros(1,cidx-1), -1, zeros(1,NR-cidx)];
NumLevels = numel(unique(JOB.model(cidx)));
if NumLevels==2 && ~isfield(JOB,'grouplabel')
JOB.grouplabel={['Non-',varnames{cidx}],varnames{cidx}};
end
M = double(JOB.model);
fname_base = varnames{cidx};
myfsl_designmatrix(fname_base, M, NumLevels, JOB.grouplabel, JOB.contrasts);

% 2-4. create mask with a percentile threshold!
unix(['ln -sf ',fname_template,' ',dir_output,'/',meastype,'_template.nii.gz']);
fname_template = [dir_output,'/',meastype,'_template.nii.gz'];
if ~isfield(JOB,'percthres')
percthres = 0.75;
else
percthres = JOB.percthres;
end
prcstr =num2str(round(percthres*100));
[a,b,c] = fileparts_gz(fname_template);
fname_mask = [a,'/',b,'_',prcstr,'perc',c];
if ~exist(fname_mask,'file')
nii = load_uns_nii(fname_template);
thrs = prctile(nii.img(~~nii.img),percthres*100);
disp(['> Mask thresholding at ',prcstr,' percentile: ',num2str(thrs)]);
nii2 = nii;
nii2.img = nii2.img > thrs;
save_untouch_nii(nii2, fname_mask);
end
unix(['slices ',fname_mask,' -o ',a,'/',b,'_',prcstr,'perc.gif']);

%% 3. just merge
fname_merged =  [dir_output,'/',meastype,'_all',fwhmstr,'.nii.gz'];
if ~exist(fname_merged,'file')
fnames_input='';
for j=1:numel(fnames)
fnames_input=[fnames_input,' ',fnames{j}];
end
unix(['fslmerge -t ',fname_merged,' ',fnames_input]);
end

%% 4. now run randomization
if ~isfield(JOB,'NumPerm')
NumPerm=10000;
else
NumPerm=JOB.NumPerm;
end
[~,qid]=myunix(['cd ',dir_glm,'; randomise_parallel -i ',fname_merged, ...
' -o ',dir_glm,fname_base,' -m ',fname_mask, ...
' -d ',dir_glm,fname_base,'.mat -t ',dir_glm,fname_base,'.con ', ...
' -n ',num2str(NumPerm),' -T --uncorrp -x ']);
qid = qid(end-4:end-1);
disp(['# condor job submitted qid=',qid]);
myunix(['waitForCONDORJobs.sh 300 ',qid]);


cd(dir0);
end
