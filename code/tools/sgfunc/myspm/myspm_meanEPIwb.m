function JOB = myspm_meanEPIwb (JOB)
if strfind(num2str(addsess), subjid)
fname_t1w=myls([dir_local,subjid,'/7T/*.SEPT/S*mp2rage*whole*UNI*_Te2.38*.nii']);
else
res=ls([dir_local,subjid,'/7T/t1w_brain.nii']);
fname_t1w = res(1:end-1)
end

trgEPI  =  myls([dir_local,subjid,'/7T/*.SEPT/S*EPI*whole*.nii']);
[p1,f1,e1]=fileparts(trgEPI);
trgEPI2 = [p1,'/m',f1,e1];
unix(['fslmaths ',trgEPI,' -Tmean ',trgEPI2]);
trgEPI = trgEPI2;
shortmag = myls([dir_local,subjid,'/7T/*.SEPT/S*field*whole*Delta*Te6*.nii']);

res=ls([dir_local,subjid,'/7T/*.SEPT/S*field*whole*Delta*Te7*.nii']);
idx = findstr(res,sprintf('\n'));
phasedif = res(idx(1)+1:end-1)

TotalReadoutTime = prepare_vdm (shortmag, phasedif, [6 7.02], trgEPI, ...
65.9196, [], fname_t1w)
correctedWholebrainEPI=[dir_local,subjid,'/7T/epi-wb.nii'];

[p1,f1,e1]=fileparts(trgEPI);
utrgEPI=[p1,'/u',f1,e1];
if strfind(num2str(addsess), subjid)
fname_t1w0p7 = [dir_local,subjid,'/7T/t1w_brain.nii'];
job1=struct('prefix','c', 'interp',1, 'name_fixed', fname_t1w0p7, ...
'name_moving', fname_t1w, 'name_others',utrgEPI);
myspm_coreg(job1);
cutrgEPI=[p1,'/cu',f1,e1];
copyfile(cutrgEPI, correctedWholebrainEPI);
else
job1=struct('prefix','c', 'interp',1, 'name_fixed', fname_t1w, ...
'name_moving', utrgEPI);
myspm_coreg(job1);
cutrgEPI=[p1,'/cu',f1,e1];
copyfile(cutrgEPI, correctedWholebrainEPI);
end

cd ([dir_local,subjid,'/7T/']);
unix('slices epi-wb.nii t1w_brain.nii -o epi-wb_over_t1w.gif');
copyfile('epi-wb_over_t1w.gif', ...
['/scr/vatikan1/skim/Tonotopy/main/fig_coreg_EPIwb/epi_over_t1w_',subjid,'.gif']);
unix('slices t1w_brain.nii epi-wb.nii -o t1w_over_epi-wb.gif');
copyfile('t1w_over_epi-wb.gif', ...
['/scr/vatikan1/skim/Tonotopy/main/fig_coreg_EPIwb/t1w_over_epi-wb_',subjid,'.gif']);
end
