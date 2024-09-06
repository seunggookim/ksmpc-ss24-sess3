function JOB = myfsl_combineclus(JOB)
% JOB = myfsl_combineclus(JOB)
%
% (cc) 2015, sgKIM.

names  = {JOB.clus.cog_name};
unames = unique(names);
nii = load_uns_nii(JOB.fname_label_nii);
img = nii.img*0;
clus_old = JOB.clus;
clus = [];
for j=1:numel(unames)
idx = find(ismember(names, unames{j}));
clus(j).id_old = [JOB.clus(idx).id];
clus(j).numVox = sum([clus_old(idx).numVox]);
clus(j).origname = unames{j};
img(ismember(nii.img, idx)) = j;
unames{j} = string2filename(unames{j});
clus(j).filename = unames{j};
end
nii.img = img;
[path1,name1,ext1] = fileparts_gz(JOB.fname_label_nii);
JOB.fname_label_nii = fullfile(path1,[name1,'_comb',ext1]);
save_untouch_nii(nii, JOB.fname_label_nii);
save(fullfile(path1,[name1,'_comb.mat']), 'clus');

JOB.clus_old = clus_old;
JOB.clus = clus;
end


%%
function fname = string2filename(str)
fname = str;
TRG={'(',')',' ',',','/'};
for i=1:numel(TRG)
fname = strrep(fname, TRG{i}, '_');
end
end
