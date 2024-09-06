function Job = myants_antsApplyTransformsTimeseries(Job)
% Job.FnameMoving
% Job.FnameFixed
% Job.Transforms
% Job.FnameOut

%% check LD_LIBRARY_PATH
%{
  NOTE: Starting up MATLAB will put MATLAB's lib paths before ANTs lib path.
  If MATLAB has its own ITK libraries (e.g., HPC) this can mask ANTs's lib.
  So I just need to switch the order in that case:
%}
LD_LIBRARY_PATH = getenv('LD_LIBRARY_PATH');
PATH_ = strsplit(LD_LIBRARY_PATH,':');
ind = contains(upper(PATH_),'ANTS');
assert(sum(ind), 'No ANTs library in LD_LIBRARY_PATH?!');
PATH_ = [PATH_(ind), PATH_(~ind)];
LD_LIBRARY_PATH = cell2mat(strcat(PATH_,':'));
setenv('LD_LIBRARY_PATH', LD_LIBRARY_PATH(1:end-1));

%% UNPACKING
tic;
[DnOut,~,~] = myfileparts(Job.FnameOut);
DnTempsrc = fullfile(DnOut, 'myantstmp', sprintf('%16i',rand*1e17));
mkdir(DnTempsrc)
logthis('Unpacking the source file "%s" into a temporary directory "%s"...\n', Job.FnameMoving, DnTempsrc)
system(sprintf('ImageMath 4 %s/.nii TimeSeriesDisassemble %s', DnTempsrc, Job.FnameMoving));

%% DATATYPE code for antsApplyTransforms
Info = niftiinfo(Job.FnameMoving);
switch (Info.Datatype)
  case 'int16'
    datatype = 'short';
  case {'single', 'int32'}
    datatype = 'float';
  otherwise
    error('FIXME!')
end

%% Transform volume by volume
files = dir(fullfile(DnTempsrc,'*.nii'));
logthis('Now transforming each of %i images...\n', numel(files))
DnTemptrg = fullfile(DnOut, 'myantstmp', sprintf('%16i',rand*1e17));
mkdir(DnTemptrg)
counter = 0;
nfiles = numel(files);
for ifile = 1:nfiles
  FnMov = fullfile(files(ifile).folder, files(ifile).name);
  FnOut = fullfile(DnTemptrg, files(ifile).name);
  cfg = struct('FnameMoving',FnMov, 'FnameFixed',Job.FnameFixed, 'FnameOut',FnOut, 'odt', datatype, ...
    'Transforms',{Job.Transforms}, 'IsCreateFig', false);
  myants_antsApplyTransforms(cfg);
  
  currentProgress = round(ifile/nfiles*10);
  if currentProgress > counter
    fprintf('..%i0%%',currentProgress)
    counter = currentProgress;
  end
end
fprintf('\n')

Info = niftiinfo(Job.FnameMoving);
Tr_sec = Info.PixelDimensions(4);

%% MERGE into a 4D volume
% This takes a bit longer but it respects the original precision
logthis('Repacking images into a single file "%s"...\n', Job.FnameOut)
setenv('FSLOUTPUTTYPE','NIFTI');
system(sprintf('fslmerge -tr %s %s/*.nii %g', Job.FnameOut, DnTemptrg, Tr_sec));
assert(isfile(Job.FnameOut))
logthis('Transformation took %.3f sec.\n', toc)

%% Remove temp directories
rmdir(DnTempsrc, 's')
rmdir(DnTemptrg, 's')
end


% function TEST()
% cfg = struct('FnameMoving',[p1,'/ua',f1,e1], ...
%   'FnameFixed','mni_funcref.nii.gz', 'Transforms',{{fn_warp,0}}, ...
%   'FnameOut',[p1,'/ua',f1,'_mni',e1]);
% myants_antsApplyTransformsTimeseries(cfg);
% end


%   
%   info = niftiinfo(job.FnameMoving);
%   switch (info.Datatype)
%     %{
%       0  :  char
%       1  :  unsigned char
%       2  :  short
%       3  :  unsigned short
%       4  :  int
%       5  :  unsigned int
%     %}
%     case 'int16'
%       typeopt = '2';
%     otherwise
%       error('FIXME!')
%   end
%   system(['ConvertImagePixelType ',fn_out,' ',fn_out,' ',typeopt ,...
%     ' 1>/dev/null']);
%   % still conmus3 data (2880 timepoints, 3mm) takes 10 min.  
