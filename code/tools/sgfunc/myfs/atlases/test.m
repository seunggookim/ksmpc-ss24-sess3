%% TEST annotation resampling:
test_resample_hcpmmp1()

%% Resample FSAVERAGE -> FSAVERAGE3,4,5,6
for i = 3:6
  resample_hcpmmp1(['fsaverage',num2str(i)])
end