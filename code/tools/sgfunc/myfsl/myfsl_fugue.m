function JOB = myfsl_fugue(JOB)

% fsl_prepare_fieldmap <scanner> <phase_image> <magnitude_image> <out_image> <deltaTE (in ms)>
[p1,f1,e1]=myfileparts(JOB.fname_mag);
json_mag = jsondecode(fileread([p1,'/',f1,'.json']));
[p2,f2,e2]=myfileparts(JOB.fname_pha);
json_pha = jsondecode(fileread([p2,'/',f2,'.json']));



end