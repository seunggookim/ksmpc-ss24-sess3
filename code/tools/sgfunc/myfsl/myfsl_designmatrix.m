function JOB=myfsl_designmatrix(JOB, fnbase, M, NumGroups, GroupLabel, Contrast)
% ## this only works for group comparison, for correlation, use myfsl_designmat.m
%
% myfsl_designmatrix(filename_prefix, M, NumGroups, GroupLabel, Contrast)
%
% makes .mat (design matrix) and .con (contrast matrix)
% for randomise/FSL for TBSS or VBM (not FMRI analysis)
%
% input:
% fnbase    = filename base
% M         = [groupA groupB covariate1 covariate2 ...]
% NumGroups = 2 (by default)
% Contrast  = [ 1 -1 zeros(1,Numwaves-2);
%              -1  1 zeros(1,Numwaves-2)] (by default)
%
% (cc) sgKIM, 2012, (solleo@gmail.com).



PPheights = range(M);
[NumPoints, NumWaves] = size(M);
if ~exist('NumGroups','var')
NumGroups = 2;
end
% contrast
Contrast  = [ 1 -1 zeros(1,NumWaves-2);
-1  1 zeros(1,NumWaves-2)];
% demeaning continuous covariates
M = M - [zeros(NumPoints, NumGroups), ...
repmat(mean(M(:,NumGroups+1:end),1),[NumPoints,1])];
% write .mat
filename = [fnbase,'.mat'];
JOB.fname_mat = filename;
fid = fopen(filename,'w');
fprintf(fid,'/NumWaves\t%d\n', NumWaves);
fprintf(fid,'/NumPoints\t%d\n', NumPoints);
fprintf(fid,'/PPheights\t');
fclose(fid);
dlmwrite(filename, PPheights, 'delimiter','\t','-append');
fid = fopen(filename,'a');
fprintf(fid,'/Matrix\n');
fclose(fid);
dlmwrite(filename, M, 'delimiter','\t','-append','precision','%.8f');
% wirte .con
filename= [fnbase,'.con'];
JOB.fname_con = filename;
fid = fopen(filename,'w');
if ~exist('GroupLabel','var')
GroupLabel={'non-Group','Group'};
end
if NumGroups == 2
fprintf(fid,['/ContrastName1\t"',GroupLabel{2},'>',GroupLabel{1},'"\n']);
fprintf(fid,['/ContrastName1\t"',GroupLabel{2},'<',GroupLabel{1},'"\n']);
elseif NumGroups == 1
fprintf(fid,['/ContrastName1\t"',GroupLabel{2},'"\n']);
end
fprintf(fid,'/NumWaves\t%d\n', NumWaves);
fprintf(fid,'/NumContrasts\t%d\n', size(Contrast,1));
% Stephan Smith said that you can ignore PPhieghts and RequiredEffect for
% TBSS (it's only for BOLD analysis using FEAT)
fclose(fid);
fid = fopen(filename,'a');
fprintf(fid,'/Matrix\n');
fclose(fid);
dlmwrite(filename, Contrast, 'delimiter','\t','-append','precision','%.8f');

end
