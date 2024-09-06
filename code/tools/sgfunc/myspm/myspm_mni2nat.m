function coord_nat = myspm_mni2nat(coord, fname_nat, fname_y)

nii=load_untouch_nii('~/MNI152_T1_2mm_brain.nii');
nii.img=nii.img*0;
ijk=round(xyz2ijk(coord,nii));
nii.img(ijk(1),ijk(2),ijk(3))=1;
fname0=[tempname,'.nii'];
save_untouch_nii(nii,fname0);

job1=[];
job1.fname_y=fname_y;
job1.fname_ref=fname_nat;
job1.fname_input=fname0;
job1.fname_interp=0;
myspm_unnorm(job1)

[p1,f1,e1]=myfileparts(fname0);
xnii=load_untouch_nii([p1,'/x',f1,e1]);
coord_nat = ijk2xyz(find3(xnii.img),xnii);
end