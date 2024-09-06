%function myspm_snpm (JOB)
load /scr/vatikan3/APConn/mat/info17_sorted.mat
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.DesignName = '2 Groups: Two Sample T test; 1 scan per subject';
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.DesignFile = 'snpm_bch_ui_TwoSampT';
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.dir = {'/scr/vatikan3/APConn/fnirt_spm/SnPM'};
%%
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans1 = {
'/scr/vatikan3/APConn/fnirt_spm/data/rcmap_fnirted_2001.nii,1'
'/scr/vatikan3/APConn/fnirt_spm/data/rcmap_fnirted_2003.nii,1'
'/scr/vatikan3/APConn/fnirt_spm/data/rcmap_fnirted_2004.nii,1'
'/scr/vatikan3/APConn/fnirt_spm/data/rcmap_fnirted_2005.nii,1'
'/scr/vatikan3/APConn/fnirt_spm/data/rcmap_fnirted_2006.nii,1'
'/scr/vatikan3/APConn/fnirt_spm/data/rcmap_fnirted_2007.nii,1'
'/scr/vatikan3/APConn/fnirt_spm/data/rcmap_fnirted_2008.nii,1'
'/scr/vatikan3/APConn/fnirt_spm/data/rcmap_fnirted_2009.nii,1'
'/scr/vatikan3/APConn/fnirt_spm/data/rcmap_fnirted_2010.nii,1'
};
%%
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans2 = {
'/scr/vatikan3/APConn/fnirt_spm/data/rcmap_fnirted_3002.nii,1'
'/scr/vatikan3/APConn/fnirt_spm/data/rcmap_fnirted_3003.nii,1'
'/scr/vatikan3/APConn/fnirt_spm/data/rcmap_fnirted_3004.nii,1'
'/scr/vatikan3/APConn/fnirt_spm/data/rcmap_fnirted_3005.nii,1'
'/scr/vatikan3/APConn/fnirt_spm/data/rcmap_fnirted_3006.nii,1'
'/scr/vatikan3/APConn/fnirt_spm/data/rcmap_fnirted_3007.nii,1'
'/scr/vatikan3/APConn/fnirt_spm/data/rcmap_fnirted_3009.nii,1'
'/scr/vatikan3/APConn/fnirt_spm/data/rcmap_fnirted_3010.nii,1'
};
%%
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(1).c = [22
20
25
24
29
25
24
38
25
26
25
23
28
28
32
24
29];
%%
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(1).cname = 'age';
%%
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(2).c = [0
0
1
1
0
1
1
0
1
1
1
0
0
1
1
0
1];
%%
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(2).cname = 'sex';
%%
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(3).c = [false
false
false
false
false
false
false
false
false
true
false
false
false
true
true
true
true];
%%
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(3).cname = 'ethnicity';
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.nPerm = 10000;
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.vFWHM = [0 0 0];
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.bVolm = 1;
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.ST.ST_later = -1;
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.masking.tm.tm_none = 1;
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.masking.im = 1;
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.masking.em = {'/scr/vatikan3/APConn/fnirt_spm/data/mean_cmap_0p1.nii,1'};
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalc.g_omit = 1;
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm.glonorm = 1;
matlabbatch{2}.spm.tools.snpm.cp.snpmcfg(1) = cfg_dep('2 Groups: Two Sample T test; 1 scan per subject: SnPMcfg.mat configuration file', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','SnPMcfg'));
matlabbatch{3}.spm.tools.snpm.inference.SnPMmat(1) = cfg_dep('Compute: SnPM.mat results file', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','SnPM'));
matlabbatch{3}.spm.tools.snpm.inference.Thr.Clus.ClusSize.CFth = 0.001;
matlabbatch{3}.spm.tools.snpm.inference.Thr.Clus.ClusSize.ClusSig.FWEthC = 0.05;
matlabbatch{3}.spm.tools.snpm.inference.Tsign = 1;
matlabbatch{3}.spm.tools.snpm.inference.WriteFiltImg.WF_no = 0;
matlabbatch{3}.spm.tools.snpm.inference.Report = 'MIPtable';
matlabbatch{4}.spm.tools.snpm.inference.SnPMmat(1) = cfg_dep('Compute: SnPM.mat results file', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','SnPM'));
matlabbatch{4}.spm.tools.snpm.inference.Thr.Clus.ClusSize.CFth = 0.001;
matlabbatch{4}.spm.tools.snpm.inference.Thr.Clus.ClusSize.ClusSig.FWEthC = 0.05;
matlabbatch{4}.spm.tools.snpm.inference.Tsign = -1;
matlabbatch{4}.spm.tools.snpm.inference.WriteFiltImg.WF_no = 0;
matlabbatch{4}.spm.tools.snpm.inference.Report = 'MIPtable';


spm_jobman('initcfg')
spm_jobman('run', matlabbatch)
