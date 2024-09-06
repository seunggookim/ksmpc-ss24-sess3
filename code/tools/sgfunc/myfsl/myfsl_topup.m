function myfsl_topup(Job)
% myfsl_topup(Job)
%
% note: "dummy scans" with interslice intensity differences are fine for
% topup
%
% Job require:
%  .fnames_pe_pair   : a pair of files with opposing phase-encoding
%                      directions 
% .fname_to_unwarp   : a file to unwarp
% (.pe_dir)          : {1x2} can be read from json for Siemens (not GE!)
% (.TotalReadoutTime): can be read from json 
% (.fname_t1w_coreg) : for visualization, coregistered T1w image
%
% Example: 
% > Scenario 1: you have 250 GE-EPI volumes in +j PE-dir (main scan) and 5
% GE-EPI volumes in -j PE-dir (additional scan) from a recent Siemens
% scanner. You have also .json files in the same directory extracted using
% dcm2niix with the same filenames as .nii files. You want to unwarp your
% main scan (250 volumes).
%
% < Solution 1: 
% myfsl_topup( strcut('fnames_pe_pair', ...
%   {'myRawEPI_+j_250volumes.nii','myRawEPI_-j_5volumes.nii'}, ...
%   'fname_to_unwarp','myRawEPI_+j_250volumes.nii') )
%
%
% > Scenario 2: you have the main scan (250 GE-EPI) in +j and one
% additional SE-EPI pair in opposing PEs (5 SP-EPI in +j; 5 SP-EPI in -j).
%
% < Solution 2: 
% myfsl_topup( strcut('fnames_pe_pair', ...
%   {'myRawSEEPI_+j_5volumes.nii','myRawSEEPI_-j_5volumes.nii'}, ...
%   'fname_to_unwarp','myRawGEEPI_+j_250volumes.nii') ) ??
%
% (cc) 2018-2022, dr.seunggoo.kim@gmail.com

%HISTORY
% 20220517: NIFTI toolbox (load_untouch_nii,...) => MATLAB (2016+) builtin
% (niftireadgz, niftiinfogz), inputs can be *.nii.gz

%% PREPARE INPUT DATA
FSLOUTPUTTYPE0 = getenv('FSLOUTPUTTYPE');
Resetter = onCleanup(@() setenv('FSLOUTPUTTYPE',FSLOUTPUTTYPE0));
setenv('FSLOUTPUTTYPE','NIFTI')

prefix_input = cell(1,2);
for i = 1:2
  [p1,f1,e1] = myfileparts(Job.fnames_pe_pair{i});
  prefix_input{i} = f1;
end
prefix_out = [prefix_input{1},'_topup'];
Job.prefix_out = [prefix_out];
fn_data = [p1,'/',prefix_out,'_data.nii'];
if ~exist(fn_data,'file')
  [img1,info1] = niftireadgz([p1,'/',prefix_input{1},e1]);
  [img2,info2] = niftireadgz([p1,'/',prefix_input{2},e1]);
  
  % take min (5 , all images) volumes:
  n = min([5 info2.ImageSize(4)]);
  
  % combine two images:
  img = cat(4, img1(:,:,:,1:n), img2(:,:,:,1:n));
  info = info1;
  info.ImageSize(4) = n*2;
  niftiwrite(img, fn_data(1:end-4), info)
  clear img1 img2 info1 info2
else
  info = niftiinfo(fn_data);
  n = info.ImageSize(4)/2;
end

%% CREATE a acquisition parameter file
fn_acq = [p1,'/',prefix_out,'_acq.txt'];
Job.fn_acq = fn_acq;

% READ json files:
json = cell(1,2);
for i=1:2
  json{i} = jsondecode(fileread(fullfile(p1,[prefix_input{i},'.json'])));
end

% DEFINE mtx:
if ~isfield(Job,'pe_dir')
  try
    mtx = cell(1,2);
    for i=1:2
      switch json{i}.PhaseEncodingDirection(1)
        % PhaseEncodingDirection: {'i','i-','j','j-','k','k-'}
        % << signs are only available for Siemens scanners
        case 'i' 
          mtx{i} = [1 0 0];
        case 'j'
          mtx{i} = [0 1 0];
        case 'k'
          mtx{i} = [0 0 1];
      end
      if numel(json{i}.PhaseEncodingDirection)==2 ...
          && strcmp(json{i}.PhaseEncodingDirection(2),'-')
        mtx{i} = -mtx{i};
      end
    end
  catch
    % PhaseEncodingAxis
    error('For GE scanners, enter .pe_dir MANUALLY!!!')
  end
else
  mtx = Job.pe_dir;
end

% WRITE a file if needed:
if ~exist(fn_acq, 'file')
  M=[
    repmat([mtx{1} json{1}.TotalReadoutTime],[n 1])
    repmat([mtx{2} json{2}.TotalReadoutTime],[n 1])
    ];
  dlmwrite(fn_acq, M, '\t')
end


%% RUN TOPUP: estimate deformation
fn_field = [p1,'/',prefix_out,'_field.nii'];
fn_unwarped = [p1,'/',prefix_out,'_unwarped.nii'];
if ~exist(fn_field,'file')
  unix(['topup --imain=',fn_data,' --datain=',fn_acq,...
    ' --out=',p1,'/',prefix_out,...
    ' --fout=',fn_field,' --iout=',fn_unwarped,' -v --fwhm=8'])
  ls(fn_field)
end


%% APPLY deformation
fn_data = Job.fname_to_unwarp;
[p1,f1,~] = myfileparts(fn_data);
fn_unwarped = [p1,'/',f1,'_unwarped.nii'];
Job.fn_unwarped = fn_unwarped;
if ~exist(fn_unwarped,'file')
  unix(['applytopup --imain=',fn_data,' --datain=',fn_acq,' --inindex=1 ',...
    '--topup=',[p1,'/',prefix_out],' --out=',fn_unwarped,' --method=jac']);
  ls(fn_unwarped)
end

fn_json = strrep(fn_data,'.nii','.json');
if isfile(fn_json)
  fn_unwarped_json = strrep(fn_json,'.json','_unwarped.json');
  copyfile(fn_json, fn_unwarped_json)
  ls(fn_unwarped_json)
end


%% GIF animation visualization
% Is the T1w coregistered to the EPI? Not really...
compare_unwarped(Job.fname_to_unwarp, fn_unwarped, Job.fname_t1w_coreg)

end
