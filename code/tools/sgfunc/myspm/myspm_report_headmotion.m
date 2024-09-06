function JOB = myspm_report_headmotion(JOB)

num_run = numel(fnames);
ax = axeslayout([num_subj num_run],'tight','tight');
for iSubj = 1:num_subj
for iRun = 1:num_run
  cd(JOB.dir_base)
  fnames= dir('rp*.txt');
  rp = load(fnames(iRun).name);
  
  axespos(ax, iRun+(iSubj-1)*num_run)
  x1=[1:size(rp,1)];
  plot(x1, [rp(:,1:3), rp(:,5:7)/pi*180])
  xlim([1 size(rp,1]);
  
end


legend({'x [mm]','y [mm]','z [mm]','r [deg]','p [deg]','y [deg]'},'location','EastOutside')
