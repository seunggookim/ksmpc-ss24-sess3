function prm = fsfast_readprm(fname_prm)
% '/scr/vatikan3/Tonotopy/FS7Tm1/2001/boldloc/001/18tones.prm'
fid=fopen(fname_prm);
C=textscan(fid,'%f %f %f %f %s');
prm.onset_sec = C{1};
prm.code = C{2};
prm.duration_sec = C{3};
prm.weight = C{4};
prm.txt = C{5};
end
