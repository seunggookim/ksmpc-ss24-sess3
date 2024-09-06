function test_resample_hcpmmp1()
% UNIT TEST of resample_hcpmmp1
% Seung-Goo KIM, 2020-05-07

fnamesTrg = resample_hcpmmp1('fsaverage');
Sides = {'lh','rh'};

for i = 1:2
  FnSrc = fullfile(pwd, 'fsaverage', [Sides{i},'.HCP-MMP1.annot']);
  [Verts0Src, LabelSrc, CotSrc] = read_annotation(FnSrc, 0);
  [Verts0Trg, LabelTrg, CotTrg] = read_annotation(fnamesTrg{i}, 0);
  
  assert(isequal(Verts0Src, Verts0Trg), 'Vertex indices mismatched')
  assert(isequal(LabelSrc, LabelTrg), 'Labels mismatched')
  assert(isequal(CotSrc, CotTrg), 'Tables mismatched')

end

fprintf('[%s:%s] ALL GOOD!\n', mfilename, datestr(now,31));

end