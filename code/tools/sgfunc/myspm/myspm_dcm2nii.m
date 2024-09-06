function myspm_dcm2nii(JOB)
% myspm_dcm2nii(JOB)
%
% JOB requires
% .fnames
% 
if ~exist('JOB','var')
 JOB=[];
end
if ~isfield(JOB,'fnames')
 [~,fnames]=myls('*');
 fnames(isempty(fnames))=[];
else
 fnames=JOB.fnames;
end
if ~isfield(JOB,'dir_out')
 [dir_out,~,~]=myfileparts(fnames{1});
else
 dir_out=JOB.dir_out;
end

matlabbatch={};
matlabbatch{1}.spm.util.import.dicom.data = reshape(fnames,[],1);
matlabbatch{1}.spm.util.import.dicom.root = 'flat';
matlabbatch{1}.spm.util.import.dicom.outdir = {dir_out};
matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;
spm('Defaults','fmri')
spm_jobman('run', matlabbatch)

%%
[~,fnames_nii]=myls([dir_out,'/*.nii']);
T=numel(fnames_nii);
if T>1
 nii0=load_untouch_nii(fnames_nii{1});
 nii0.img=repmat(nii0.img,[1 1 1 T]);
 nii0.hdr.dime.dim(1)=4;
 nii0.hdr.dime.dim(5)=T;
 for t=2:T
  nii=load_untouch_nii(fnames_nii{t});
  nii0.img(:,:,:,t) = nii.img;
 end
 hdr=spm_dicom_headers(fnames{1});
 nii0.hdr.dime.pixdim(5)=hdr{1}.RepetitionTime/1000;
end
[~,f1,~]=myfileparts(dir_out);
save_untouch_nii(nii0,[dir_out,'/',f1,'.nii']);
for t=1:T
 delete(fnames_nii{t})
end

end