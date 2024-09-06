function JOB = myfsl_prbtrk_bidir (JOB)

[~,name_roi1,~] = myfileparts(JOB.fname_roi1);
[~,name_roi2,~] = myfileparts(JOB.fname_roi2);

[p1,f1,e1] = myfileparts(JOB.dname_dwi);
dn_probtrackx = [p1,'/',f1,e1,'.prbtrk.bidir.',JOB.name_probtrackx];

job1=[];
job1.dname_dwi = JOB.dname_dwi;
job1.dname_probtrackx = [dn_probtrackx,'/',name_roi1,'-to-',name_roi2];
job1.fname_seed = JOB.fname_roi1;
job1.fname_waypoint = JOB.fname_roi2;
job1.fname_stop = JOB.fname_roi2;
job1.optarg = ' --usef=0.2';
job1 = myfsl_probtrackx2 (job1);

job2=[];
job2.dname_dwi = JOB.dname_dwi;
job2.dname_probtrackx = [dn_probtrackx,'/',name_roi2,'-to-',name_roi1];
job2.fname_seed = JOB.fname_roi2;
job2.fname_waypoint = JOB.fname_roi1;
job2.fname_stop = JOB.fname_roi1;
job2.optarg = ' --usef=0.2';
job2 = myfsl_probtrackx2 (job2);

numtracks1 = load_untouch_nii(job1.fname_numtracks);
lengths1 = load_untouch_nii(job1.fname_lengths);

numtracks2 = load_untouch_nii(job2.fname_numtracks);
lengths2 = load_untouch_nii(job2.fname_lengths);

cd(dn_probtrackx)
unix(['echo ',JOB.fname_roi1,' > roilist.txt']);
unix(['echo ',JOB.fname_roi2,' >> roilist.txt']);

unix(['cat ',job1.dname_probtrackx,'/waytotal > waytotal']);
unix(['cat ',job2.dname_probtrackx,'/waytotal >> waytotal']);

numtracks = numtracks1;
numtracks.img = (numtracks1.img + numtracks2.img);
save_untouch_nii(numtracks,'fdt_paths.nii.gz')

brain = load_untouch_nii([JOB.dname_dwi,'/bmnodiff_unwarped.nii']);

cfg=struct('contournum',3, 'contourcolor',[.5 .5 .5],...
  'colormap',hot, 'caxisprc',[0 100]);
cfg.colorbartitle='# tracks [smp]';
cfg.fname_png = 'numTracks.png';
slices(numtracks.img, cfg, brain.img);

lengths = lengths1;
lengths.img = 0.5*(lengths1.img + lengths2.img);
save_untouch_nii(lengths,'fdt_paths_lengths.nii.gz')

cfg.colorbartitle='mean dist [mm]';
cfg.colormap=[0 0 0; parula];
cfg.fname_png = 'meanDist.png';
slices(lengths.img, cfg, brain.img);

numtracksbylengths = numtracks;
numtracksbylengths.img = numtracks1.img .* lengths1.img ...
  + numtracks2.img .* lengths2.img;
save_untouch_nii(numtracksbylengths,'fdt_paths_by_lengths.nii.gz')

cfg.colorbartitle='# tracks x mean dist [smp.mm]';
cfg.colormap=hot;
cfg.fname_png = 'numTracksByDist.png';
slices(numtracksbylengths.img, cfg, brain.img);

end

%%
