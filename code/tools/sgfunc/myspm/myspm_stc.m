function myspm_stc(job)
% job.fname_epi
% job.NumFrames
% job.TR_sec
% job.slice_order_msec
% job.ref_slice_msec
%
% separated for parallel jobs


st1.scans={};
for t = 1:job.NumFrames
  st1.scans{1}{t,1} = [job.fname_epi,',',num2str(t)];
end
st1.nslices = numel(job.slice_order_msec);
st1.tr = job.TR_sec;
st1.ta = job.TR_sec-(job.TR_sec/st1.nslices); % NOTE: this value will
% not be used because .so and .refslice are in msec (see SPM Batch
% Editor help)
st1.so = job.slice_order_msec;
st1.refslice = job.ref_slice_msec;
st1.prefix = 'a';
matlabbatch = {};
matlabbatch{1}.spm.temporal.st = st1;
[p1,f1,~] = myfileparts(job.fname_epi);
fname_matlabbatch = [p1,'/myspm_fmriprep_stc_',f1,'.mat'];
save(fname_matlabbatch,'matlabbatch');
spm_jobman('initcfg')
spm_jobman('run', matlabbatch);
end
