function JOB = myfsl_probtrkx (JOB)
% JOB = myfsl_probtrkx (JOB)
%
% JOB requires:
%  .coords_mm  {1xN} for N subjects, contaitning <Px3> coordinates for P voxels
% or
%  .fname_seed '1xN' seed
%
%  .fname_wpt  {1xN} way-point
%  .fname_trg  {1xN} target
%  .fname_stp  {1xN} stop-point
%
% (cc) 2015-2019, sgKIM. solleo@gmail.com  https://ggooo.wordpress.com


dir_base = JOB.dir_base;
dir_old = pwd;
if ~isfield(JOB,'dir_sub'),    JOB.dir_sub='fdt_paths'; end
if ~isfield(JOB,'dir_tbss'),   JOB.dir_tbss=dir_base;   end
if ~isfield(JOB,'fname_seed') && ~isfield(JOB,'coords_mm') ...
    && ~isfield(JOB,'fname_label_mat')
  error('Enter either .coords_mm or .fname_label_mat!');
end
issimple = isfield(JOB,'coords_mm');
islabel  = isfield(JOB,'fname_label_mat');

if islabel
  x = load(fullfile(JOB.fname_label_mat),'clus'); % why not? you have seed indices!
  JOB.clus = x.clus;
  for c=1:numel(JOB.clus)
    JOB.clus(c).filename=str2filename(JOB.clus(c).cog_name);
  end
end

subjid = fsss_subjID(JOB.subjID);
subjid = subjid{1};
dir1   = fullfile(dir_base,subjid);
dir_bedpostx = [dir1,'.bedpostX/'];
dir_probtrkx = [dir1,'.probtrkX/',JOB.dir_sub];
[~,~]=mkdir(dir_probtrkx);
% regarg  = fsss_fs2native(struct('subjid',subjid, 'meastype', modal, ...
%   'fname_vol',fname_like, 'fsdir',dir_fs),1);
prefix_samples = [dir_bedpostx,'merged'];
fname_mask     = [dir_bedpostx,'nodif_brain_mask.nii.gz'];

% coordinates
trkmode='';
if issimple
  fname_seed = fullfile(dir_probtrkx,JOB.name_seed); % in voxel index!!
  coords_mm = JOB.coords_mm{i};
  M = coords_mm*0;
  for j=1:size(coords_mm,1);
    M(j,:)=xyz2ijk(coords_mm(j,:), fname_mask);
  end
  rM = round(M);
  dlmwrite(fname_seed,rM,'delimiter','\t');
  trkmode=' --simple ';
  
  % clusters from TBSS (vmap_???.nii for each mask file)
elseif islabel
  JOB.fname_label_nii = [JOB.dir_tbss,'/FA/',subjid,'_FA_clusidmap.nii.gz'];
  nii = load_uns_nii(JOB.fname_label_nii);
  numseeds = numel(JOB.clus);
  FNAME_SEEDS=cell(1,numseeds);
  for j=1:numseeds
    seed=nii;
    seed.img = double(round(nii.img) == j);
    fname_seed = fullfile(dir_probtrkx, ...
      ['seed',pad(j,3),'_',JOB.clus(j).filename,'.nii.gz']);
    seed.hdr.dime.datatype=2;
    save_untouch_nii(seed, fname_seed);
    unix(['mri_convert --like ',fname_mask,' ',fname_seed,' ',fname_seed]);
    unix(['fslmaths ',fname_seed,' -thr 0.01 -bin ',fname_seed]);
    FNAME_SEEDS{j} = fname_seed;
    nii_seed_bin = load_uns_nii(fname_seed);
    if j==1, nii_nn = nii_seed_bin; end
    nii_nn.img(~~nii_seed_bin.img) = j;
  end
  save_untouch_nii(nii_nn, fullfile(dir_probtrkx, 'label_nn.nii.gz'));
  trkmode='';
  
  % or binary masks
else % single mask
  if ~iscell(JOB.fname_seed)
    [~,~,ext1]=fileparts_gz(JOB.fname_seed);
    %     if sum(ismember(ext1,{'nii.gz','nii','img','hdr'}))
    % for a seed mask
    if ~strcmp('/',JOB.fname_seed(1))
      src = fullfile(dir1,JOB.fname_seed);
    else
      src = JOB.fname_seed;
    end
    [~,~,ext1] = fileparts_gz(JOB.fname_seed);
    trg = fullfile(dir_probtrkx,['seed',ext1]);
    unix(['ln -s ',src,' ',trg]);
    fname_seed = trg;
    trkmode='';
    
  else % multiple masks
    % create seed list
    fid = fopen([dir_probtrkx,'/seedlist.txt'], 'w');
    for k=1:numel(JOB.fname_seed)
      % find absolute path
      if ~strcmp('/',JOB.fname_seed{k}(1))
        fname_seed = fullfile(dir1,JOB.fname_seed{k});
      else
        fname_seed = JOB.fname_seed{k};
      end
      fprintf(fid,'%s\n', fname_seed);
    end
    fclose(fid);
    fname_seed = [dir_probtrkx,'/seedlist.txt'];
    trkmode=' --network  ';
  end
end

% if a target mask is given
if isfield(JOB,'fname_trg')
  % create 'target mask list' file
  fid = fopen([dir_probtrkx,'/trglist.txt'], 'w');
  for k=1:numel(JOB.fname_trg)
    % find absolute path
    if ~strcmp('/',JOB.fname_trg{k}(1))
      fname_trg = fullfile(dir1,JOB.fname_trg{k});
    else
      fname_trg = JOB.fname_trg{k};
    end
    fprintf(fid,'%s\n', fname_trg);
  end
  fclose(fid);
  trgarg = [' --targetmasks=',dir_probtrkx,'/trglist.txt '];
else
  trgarg = '';
end

% if a stop mask is given
if isfield(JOB,'fname_stp')
  if iscell(JOB.fname_stp)
    fname_stp=[dir_probtrkx,'/stopmask.nii.gz'];
    cmd=['FSLOUTPUTTYPE=NIFTI_GZ; fslmaths '];
    for j=1:numel(JOB.fname_stp)
      cmd=[cmd ,JOB.fname_stp{j}, ' -add '];
    end
    cmd(end-5:end)=[];
    cmd=[cmd ' -bin ',fname_stp];
    unix(cmd);
  end
  stparg = [' --stop=',fname_stp];
else
  stparg = '';
end

% if a waypoint mask (get only streamline that touch any of waypoints) is given
if isfield(JOB,'fname_wpt')
  % create 'waypoint mask list' file
  fid = fopen([dir_probtrkx,'/wptlist.txt'], 'w');
  for k=1:numel(JOB.fname_wpt)
    % find absolute path
    if ~strcmp('/',JOB.fname_wpt{k}(1))
      fname_wpt = fullfile(dir1,JOB.fname_wpt{k});
    else
      fname_wpt = JOB.fname_wpt{k};
    end
    fprintf(fid,'%s\n', fname_wpt);
  end
  fclose(fid);
  wptarg = [' --waypoints=',dir_probtrkx,'/wptlist.txt --waycond=AND'];
else
  wptarg = '';
end

setenv('FSLOUTPUTTYPE','NIFTI_GZ')
cd(dir_probtrkx);

cmd=['probtrackx2 ',trkmode,' --forcedir -V 1 ',...
  ' -l -c 0.2 -S 2000 --steplength=0.5 --fibthresh=0.01 ', ...
  ' --distthresh=0.0 --sampvox=0.0 ', ...
  ' -s ',prefix_samples,' -m ',fname_mask,' -x ',fname_seed,...
  ' -o vmap --dir=',dir_probtrkx,' ',...
  ' --opd -P 10000 --rseed=',num2str(randi(99999)), ...
  ' ',trgarg,' ',wptarg,' ',stparg];
unix(cmd);

% now need to rename results and move to subdirectory
if issimple
  for j=1:size(rM)
    src=sprintf(['%s/vmap_%i_%i_%i.nii.gz'],dir_probtrkx,rM(j,:));
    trg=sprintf(['%s/vmap_%s_%i_%i_%i.nii.gz'],dir_probtrkx,JOB.SEED_NAMES{j},rM(j,:));
    movefile(src,trg);
  end
elseif islabel
  movefile('vmap.nii.gz', ['vmap_',JOB.clus(j).filename,'.nii.gz'])
end

cd(dir_old);
end
