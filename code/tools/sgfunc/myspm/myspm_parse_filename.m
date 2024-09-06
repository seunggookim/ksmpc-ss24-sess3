function [fwhm,isW,isU,isA,origname] = myspm_parse_filename(fname)
% [fwhm,isW,isU,isA,origname] = myspm_parse_filename(fname)
%
% only works for a speicific preproc order but allows some omissions
% STC >> unwarp >> norm >> smoothing
%
% (cc) 2018, sgKIM.

[~,fname,~]=myfileparts(fname);
fwhm=nan;
if strcmp(fname(1),'s')
  fwhm=0;
  fname(1)=[];
  n=[];
  while isnumeric(str2double(fname(1))) && ~isnan(str2double(fname(1)))
    n=[n fname(1)];
    fname(1)=[];
  end
  fwhm=str2double(n);
end
if strcmp(fname(1),'w')
  isW=true;
  fname(1)=[];
else
  isW=false;
end
if strcmp(fname(1),'u')
  isU=true;
  fname(1)=[];
else
  isU=false;
end
if strcmp(fname(1),'a')
  isA=true;
  fname(1)=[];
else
  isA=false;
end
origname=fname;

end