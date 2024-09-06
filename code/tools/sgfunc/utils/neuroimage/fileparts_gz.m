function [path,name,ext]=fileparts_gz(filename)
% [path,name,ext]=fileparts_gz(filename)
% (cc) 2014, sgKIM.

[path,name,ext] = fileparts(filename);
if strcmp(ext,'.gz')
ext = [name(end-3:end),ext];
name = name(1:end-4);
end

end
