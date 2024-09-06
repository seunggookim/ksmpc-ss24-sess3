function JOB = myspm_rsfc_seed_vox2vox (JOB)
% JOB = myspm_rsfc_seed_vox2vox (JOB)
%
% computes correlation matrix between voxels in region A vs. voxels in
% region B
%
% JOB requires:
%  .SEEDPAIRS;
%  .dir_base
%
% (cc) 2015, sgKIM   mailto://solleo@gmail.com
global overwrite
if isempty(overwrite), overwrite=0; end
if ~isfield(JOB,'aparcname'), JOB.aparcname='aparc.a2009s+aseg'; end
dir0=pwd;
if ~isfield(JOB,'fwhm'), JOB.fwhm=0; end
fwhm = sprintf('%0.1f',JOB.fwhm);
%%
JOB=[];
JOB.dir_base='/scr/vatikan3/APConn/rest12.410/';
JOB.aparcname='aparc+aseg';
side={'lh-','rh-'};
JOB.SEEDPAIRS={  {'transversetemporal','superiortemporal'} };
JOB.subjID = 2001;
JOB.fwhm=4.6;
fwhm = sprintf('%0.1f',JOB.fwhm);

subjID = fsss_subjID(JOB.subjID);
for j=1:numel(JOB.SEEDPAIRS)
%seedpairs = JOB.SEEDPAIRS{j};
%   dir_mni=[JOB.dir_base,'/mni152/s',fwhm,'_native_mean_',seedname];
%   [~,~]=mkdir(dir_mni);
%%
for n=1:numel(subjID)
subjid=subjID{n};
%     if ~exist([dir_mni,'/wcbeta2_',subjid,'.nii'],'file')||overwrite
dir1=[JOB.dir_base,'/',subjid];
cd (dir1);

% smoothing
job1=[];
job1.fnames=[dir1,'/fruarest410.nii'];
job1.fwhm = JOB.fwhm;
sfname = [dir1,'/s',fwhm,'fruarest410.nii'];
if ~exist(sfname,'file');
myspm_smooth(job1);
end
P = spm_vol(sfname);
y = load_uns_nii(sfname,1); % first frame

hf=figure('position',[1064         605         730         535]);
for s=1:2
seedpairs = { [side{s},JOB.SEEDPAIRS{j}{1}], [side{s},JOB.SEEDPAIRS{j}{2}]};
% Find seed coordinates (IJK)
aparc = load_uns_nii([dir1,'/o',JOB.aparcname,'.nii']);
[~,~,roiname,roi_label] ...
= load_roiname([dir1,'/o',JOB.aparcname,'1.txt']);

for a=1:2
label{a} = roi_label(ismember(roiname, seedpairs{a}));
idx{a} = ismember(aparc.img, label{a});
IJK{a} = find3(idx{a});
Y{a} = spm_get_data(P, IJK{a}');
end
x = pca([IJK{1}; IJK{2}],1); % it should be from ventral to dorsal

paramaxis{1} = x(1:size(IJK{1},1));
paramaxis{2} = x(size(IJK{1},1)+1:end);
C = corr(Y{1},Y{2});

subplot(2,1,s)
hold on;
vecx=[]; vecy=[];
for c=1:size(Y{1},2)
%scatter(paramaxis{2}-paramaxis{1}(c),C(c,:),'k.')
vecx=[vecx; paramaxis{2}-paramaxis{1}(c)];
vecy=[vecy;C(c,:)'];
end
jhist(vecx,vecy); colormap(flipud(gray)); ylim([-1 1])
xlabel('dist between voxels along the PC#1 of STG coords (mm)'); ylabel('corr');
title([subjid,':',seedpairs{1},' vs. ',seedpairs{2}])
end
screen2png(['/scr/vatikan3/APConn/figures.local/PAC_vs_STG_corr/corr_over_dist_',subjid,'.png']);
close(hf);

% and.. and classify????
end


%%
end



end


