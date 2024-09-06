function minp = myfsl_minp(fname_tfce)
[~,minp] = unix(['fslstats ',fname_tfce,' -R']);
minp = textscan(minp,'%f %f');
minp = 1 - (minp{2});
end
