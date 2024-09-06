function JOB = myspm_long12_pairwise (JOB)
% JOB = myspm_long12_pairwise (JOB)
%
% JOB requires:
% .fname1
% .fname2
%(.vox_mm)
%
% (cc) 2017, sgKIM, solleo@gmail.com

%% 1. Creating a subject-template (and JD from TP1 to TP2)
disp('[1] Creating a subject-template..');
pw=[];
pw.vols1={JOB.fname1};
pw.vols2={JOB.fname2};
pw.tdif = 1;
pw.noise = NaN;
pw.wparam = [0 0 100 25 100];
pw.bparam = 1000000;
pw.write_avg = 1;
pw.write_jac = 1;
pw.write_div = 0;
pw.write_def = 1;
matlabbatch={};
matlabbatch{1}.spm.tools.longit{1}.pairwise = pw;
[p1,f1,e1]=fileparts(pw.vols1{1});
[p2,f2,e2]=fileparts(pw.vols2{1});
if isempty(p1),p1=pwd;end
fname_out=[p1,'/avg_',f1,e1];
if ~exist(fname_out,'file')
 save([p1,'/spm12_long-pairwise.mat'], 'matlabbatch');
 spm_jobman('initcfg')
 spm_jobman('run', matlabbatch)
end
isdone(fname_out,1)

%% 2. Unified segmentation the subj-template to MNI152
disp('[2] Running unified segmentation on the subj-template..');
% Unified segmentation
fname_in = fname_out;
fname_out=[p1,'/y_avg_',f1,e1];
if ~exist(fname_out,'file')
 myspm_seg12(fname_in,1);
 myunix(['mv ',p1,'/*c?avg_',f1,e1,' /tmp/']);
 myunix(['mv ',p1,'/mavg_',f1,e1,' /tmp/']);
end
isdone(fname_out,2)

%% 3. Normalized T1w for sanity check
disp('[3] Normalizing the JD into MNI152..');
fname_out=[p1,'/wavg_',f1,e1];
if ~exist(fname_out,'file')
 job2=JOB;
 job2.fname_moving=[p1,'/avg_',f1,e1];
 job2.fname_deform=[p1,'/y_avg_',f1,e1];
 myspm_norm(job2);
end
isdone(fname_out,3)

%% 4. Registering the JD into MNI152
disp('[4] Normalizing the JD into MNI152..');
job2=JOB;
job2.fname_moving=[p1,'/jd_',f1,'_',f2,e1];
job2.fname_deform=[p1,'/y_avg_',f1,e1];
fname_out=[p1,'/wjd_',f1,'_',f2,e1];
if ~exist(fname_out,'file')
 myspm_norm(job2);
end
isdone(fname_out,4)

end

function isdone(fname_out,proci)
if exist(fname_out,'file')
 disp(['[',num2str(proci),'] done: ',fname_out])
else
 error(['[',num2str(proci),'] failed: ',fname_out,' was not created'])
end
end