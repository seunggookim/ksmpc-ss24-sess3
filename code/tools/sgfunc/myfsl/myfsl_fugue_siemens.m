function JOB = myfsl_fugue_siemens(JOB)
% JOB requires:
% .fname_mag
% .fname_pha


% fsl_prepare_fieldmap <scanner> <phase_image> <magnitude_image> <out_image> <deltaTE (in ms)>
[p1,f1,e1]=myfileparts(JOB.fname_mag);
json_mag = jsondecode(fileread([p1,'/',f1,'.json']));
[p2,f2,e2]=myfileparts(JOB.fname_pha);
json_pha = jsondecode(fileread([p2,'/',f2,'.json']));
delaTE_ms = (json_pha.EchoTime-json_mag.EchoTime)*1000
fname_mag_brain=[p1,'/bm',f1,e1];
if ~exist(fname_mag_brain,'file')
  myspm_seg12(JOB.fname_mag,'ss')
end
fname_fmap = ['fmaprad_',f1,e1];
if ~exist(fname_fmap,'file')
  unix(['fsl_prepare_fieldmap SIEMENS ',JOB.fname_pha,' ',fname_mag_brain,' ',fname_fmap,' ',num2str(delaTE_ms)])
end

% apply fieldmap to EPI using FEAT


end