function JOB=myfsl_randomise(JOB)
% examples of the fields of JOB:
% JOB.dir_name=['/scr/vatikan3/Roberta/TPIAMO/dparsf/GLMs/meanFD_s8_',ROI_NAMES{i},'_',BHV_NAMES{j}];
% JOB.fname_all=/fullpath/of/your/skeletonised/measures.nii.gz;
% JOB.vn.val, JOB.vn.name (variable of nuissance)
% JOB.vi.val, JOB.vi.name (variable of interest)
% JOB.masking = full path for an explicit mask
% JOB.num_rand = <integer> # of randomizations
%
% (cc) 2015. sgKIM. solleo@gmail.com

fname_all = JOB.fname_all;
fname_mask = JOB.fname_mask;

mkdir(JOB.dir_name)
output_prefix =[JOB.dir_name,'/',JOB.analysis];
fname_design  =[JOB.dir_name,'/design.mat'];
N = numel(JOB.vi.val);
if size(JOB.vi.val,1) ~=N
JOB.vi.val = JOB.vi.val';
end
M = [ones(N,1), JOB.vi.val];
if isfield(JOB,'vn')
K = numel(JOB.vn);
MM = zeros(N,K);
for c=1:K
MM(:,c) = JOB.vn(c).val;
end
M = [M MM];
else K=0;
end
dlmwrite(fname_design, M, '\t');
fname_header = [JOB.dir_name,'/design.header'];
fid = fopen(fname_header,'w');
fprintf(fid, '/NumWaves %i\n', K+2);
fprintf(fid, '/NumPoints %i\n', N);
fprintf(fid, '/Matrix\n');
fclose(fid);
system(['cat ',fname_header,' ',fname_design,'> /tmp/mat.txt']);
movefile('/tmp/mat.txt', fname_design);

fname_contrast=[JOB.dir_name,'/design.con'];
if isfield(JOB,'vn')
M = [0  1 zeros(1,numel(JOB.vn)); 0 -1 zeros(1,numel(JOB.vn))];
else
M = [0 1; 0 -1];
end
dlmwrite(fname_contrast, M, '\t');
fname_header = [JOB.dir_name,'/contrast.header'];
fid = fopen(fname_header,'w');
fprintf(fid, '/ContrastName1 %s\n', K+2);
fprintf(fid, '/ContrastName2 %i\n', N);
fprintf(fid, '/NumWaves %i\n', K+2);
fprintf(fid, '/NumContrasts 2\n');
fprintf(fid, '/Matrix\n');
fclose(fid);
system(['cat ',fname_header,' ',fname_contrast,'> /tmp/mat.txt']);
movefile('/tmp/mat.txt', fname_contrast);

num_rand=JOB.num_rand; % for test
system(['randomise_parallel -i ',fname_all,' -o ',output_prefix,' -m ',fname_mask,...
' -d ',fname_design,' -t ',fname_contrast,' -n ',num2str(num_rand),' --T2'])

% find clusters

% find p-values

% create a table

% create figures

end
