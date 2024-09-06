function myfs_addpath
[mypath,~,~] = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(mypath,'external')))
addpath(genpath(fullfile(mypath,'utilities')))
end
