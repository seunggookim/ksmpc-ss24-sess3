function [lastfname,allfnames] = myls(str,isDir)
% [lastfname,allfnames] = myls(str,isDir)
if ~exist('isDir','var'), isDir=''; end
try
  if isDir
    res=ls('-d',str);
  else
    res=ls(str);
  end
catch ME
  lastfname={}; allfnames={}; return    
end
idx=strfind(res,sprintf('\n'));
idx = sort([idx strfind(res,sprintf('\t'))]);
if numel(idx)>1
  lastfname = res(idx(end-1)+1:end-1);
else
  lastfname = res(1:end-1);
end

allfnames=cell(1,numel(idx));
if nargout>1
  idx = [0 idx];
  for j=2:numel(idx)
    allfnames{j} = res(idx(j-1)+1:idx(j)-1);
  end
end
allfnames(1)=[];
%ls(lastfname)
end
