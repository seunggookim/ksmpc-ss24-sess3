function myspm_mywds(fname_epi_orig)
% JOB requires:
%
% .fname_epi_orig


%{ 
Let's try Patel's way first as described in
Patel & Bullmore, 2016, "A wavelet-based estimator of the degrees of
freedom in denoised fMRI time series for probabilistic testing of
functional connectivity and brain graphs", NeuroImage.

Core Image Processing included the following steps: 
(i) slice ac- quisition correction; 
(ii) rigid-body head movement correction to the first frame of data; 
(iii) obliquity transform to the structural image; 
(iv) affine co-registration to the skull-stripped structural image using 
a gray matter mask; 
(v) standard space transform to the MNI152 template in MNI space; 
(v) spatial smoothing (6 mm full width at half maximum); 
(vi) a within-run intensity normalization to a whole-brain median of 1000. 
All spatial transforms were applied in one step to avoid incremental blurring 
of the data that can occur from multiple independent transforms.

Denoising steps included: 
(vii) wavelet despiking (performed voxel- wise with the BrainWavelet Toolbox); 
(viii) confound signal regression including the 6 motion parameters estimated 
in (ii), their first order temporal derivatives, and ventricular cerebrospinal 
fluid (CSF) signal; and 
(ix) a wavelet “band-pass” filter. In this last step, the Maximal Over- 
lap Discrete Wavelet Transform (MODWT) was used to produce a set of scales 
(frequency bands), and coefficients from scales representing fre- 
quency bands of interest were recomposed to produce frequency- 
filtered time series.

%}

[p1,f1,~]=myfileparts(fname_epi_orig);
epi_suffix=f1;
cd(p1)

% 1. (viii) regression first
fname_rp=['rp_a',epi_suffix,'.txt'];
rp=load(fname_rp);
% rp=[...
%  rp(:,1:3), l2norm([0 0 0; diff(rp(:,1:3))]) , ...
%  rp(:,4:6), l2norm([0 0 0; diff(rp(:,4:6))]) ];
% 1-2. compcor regressors (3) 
fname_epi=['ua',epi_suffix,'.nii'];
job2=[];
job2.dir_data=p1;
job2.fname_epi=fname_epi;
hdr=load_untouch_header_only([epi_suffix,'.nii']);
job2.TR_sec=hdr.dime.pixdim(5);
disp(['TR = ',num2str(job2.TR_sec),' sec']);
job2.num_pcs=1;
job2.onlycsf=1;
job2=myy_compcor(job2);
cc=load(job2.fname_cc);
% 3-3. get residual
nii=load_untouch_nii(fname_epi);
y=img2y(nii.img);
M=[ones(size(y(:,1))) rp ];
yres = y-M*(pinv(M'*M)*M'*y);
nii.img = y2img(yres, size(nii.img));
save_untouch_nii(nii,['r1',fname_epi]);

% 1. (vii) WDS

fname_wds=['Wr1',fname_epi];
if ~exist(fname_wds,'file')
 WaveletDespike(['r1',fname_epi],'WDS');
 unix(['gunzip WDS*.nii.gz']);
 unix(['mv WDS_wds.nii ',fname_wds])
end
ls(fname_wds)

% 2. need to swap matrices and change vox-to-world matrix...
job1=[];
job1.fname_moving = fname_wds;
job1.fname_fixed  = ['meanua',epi_suffix,'.nii'];
job1.interp=0;
myspm_coreg4d(job1)
fname_wds=['o',fname_wds];


% 4. 

% 5. 

end
