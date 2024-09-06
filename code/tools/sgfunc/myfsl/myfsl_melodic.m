function JOB=myfsl_melodic(JOB)

path1='/home/MINTS/sgk26/Dropbox/sgfunc/myfsl';
[path1,~,~]=fileparts(mfilename);
C=textscan(fopen(fullfile(path1,'template_melodic.fsf')),'%s\n','delimiter','');
CMD=C{1};
fid=fopen('~/test.fsf','w');
for j=1:numel(CMD)
 
 case
 
 fprintf(fid,'%s\n',CMD{j});
 
end
fclose(fid)

end