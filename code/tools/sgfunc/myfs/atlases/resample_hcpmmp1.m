function FnamesTrg = resample_hcpmmp1(TrgSubid)
% FnamesTrg = resample_hcpmmp1(TrgSubid)
% Seung-Goo KIM, 2020-05-07

SrcSubid = 'fsaverage';
% SubDir = getenv('SUBJECTS_DIR');

Sides = {'lh','rh'};
FnamesTrg = cell(1,2);
for iHemi = 1:2
  FnSrc = fullfile(pwd, SrcSubid, [Sides{iHemi},'.HCP-MMP1.annot']);  
  
  DnTrg = fullfile(pwd, TrgSubid);
  if ~isfolder(DnTrg)
    mkdir(DnTrg)
  end
  if strcmpi(SrcSubid, TrgSubid)
    warning('Identical subjects: for TESTING')
    FnTrg = fullfile(pwd, TrgSubid, [Sides{iHemi},'.HCP-MMP1.test.annot']);
  else
    FnTrg = fullfile(pwd, TrgSubid, [Sides{iHemi},'.HCP-MMP1.annot']);
  end
  
  cmd = ['mri_surf2surf --srcsubject ',SrcSubid,' --sval-annot ',FnSrc,...
    ' --trgsubject ',TrgSubid,' --tval ',FnTrg,' --hemi ',Sides{iHemi},...
    ' --mapmethod nnf'];
  system(cmd);
  
  FnamesTrg{iHemi} = FnTrg;
end

% [Verts0Trg, LabelTrg, CotTrg] = read_annotation(FnTempTrg, 0);


end
