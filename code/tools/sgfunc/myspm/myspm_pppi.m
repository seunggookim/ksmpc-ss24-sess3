function JOB = myspm_pppi (JOB)
% JOB = myspm_pppi (JOB)
%
% JOB requires:
% .subjid
% .dir_glm
% .fname_voi
% .voi_name
% .convec
% .conname
%
% (cc) sgKIM, 2017.
cd(JOB.dir_glm)
%% Check whether VOI fits into the mask.nii
fname1=mydir('mask*');
mask=load_uns_nii(fname1);
voi =load_uns_nii(JOB.fname_voi);
if sum(double(voi.img(:))) ~= sum(double(voi.img(:)).*double(mask.img(:)))
 error('The VOI mask includes voxels outside the GLM mask! Use VOI mask that fits to interaction of all GLM masks!');
end

%%
P = [];
P.subject   = JOB.subjid;
P.directory = JOB.dir_glm;
P.VOI       = JOB.fname_voi;
P.Region    = JOB.voi_name;%STR{j};
P.estimate  = 1;
P.contrast  = 0;
P.analysis  = 'psy';
P.extract   = 'eig';
P.method    = 'cond';
P.weighted = 0;
P.contrast = 0;
P.Estimate = 1;
PPPI(P)
% delete ResI files (because I modified my SPM12 to leave ResI by me.
dir_ppi = [P.directory,'/PPI_',P.Region];
myunix(['rm -f ',dir_ppi,'/ResI*.nii']);
% move all the log/mat/txt
[~,~]=mkdir([dir_ppi,'/pppi_log/']);
myunix(['mv ',P.directory,'/',P.subject,'_PP* ',dir_ppi,'/pppi_log/']);
myunix(['mv ',P.directory,'/',P.subject,'*',JOB.voi_name,'* ',dir_ppi,'/pppi_log/']);

% compute contrast files (requires .conname and .convec)
if isfield(JOB,'convec') && isfield(JOB,'conname')
JOB2=JOB;
JOB2.dir_glm=dir_ppi;
h1=figure;
myspm_mycon(JOB2)
close(h1)
end

end
