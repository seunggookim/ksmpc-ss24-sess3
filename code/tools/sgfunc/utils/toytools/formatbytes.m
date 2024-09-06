function r = formatbytes(bytes)
%formatbytes returns a human-readable formatted string of a given bytes
%
% A modified version of: http://gaidi.ca/weblog/human-readable-bytes-in-matlab
% BASE 2 is for KiB (kibi-), MiB (mebi-), GiB (gibi-), TiB (tebi-) :IEC
% 80000-13
% 
% (cc) 2021, sgKIM.

s = {'B', 'KiB', 'MiB', 'GiB', 'TiB'}; 
e = min(floor(log(bytes)/log(1024)), 4);
r = sprintf('%7.2f %s', (bytes/(1024^floor(e))), s{e+1});
end
