function tab = myspm_readtable(fname)
if ~exist('fname','var')
  fname=mydir('spm*.csv');
  if ~ischar(fname)
    fname=fname{end};
  end
end

try tab=readtable(fname,'delimiter','\t');
catch
fid = fopen(fname,'r');
tab_fmt = '%s %s %f %s %f %f %f %f %f %f %f %f %f %s %f';
C = textscan(fid,tab_fmt ,'headerlines',1, 'delimiter','\t');
fclose(fid);

tab=[];
tab.effectSize=C{3};
tab.stats=C{4};
tab.peakT=C{5};
tab.peakZ=C{6};
tab.uncor_pval=C{7};
tab.cor_peak_pval=C{8};
tab.cor_clus_pval=C{9};
tab.extent_voxel=C{10};
tab.mni_xyz=[C{11} C{12} C{13}];
tab.struct_name=C{14};
tab.struct_name_prob=C{15};
end
end
