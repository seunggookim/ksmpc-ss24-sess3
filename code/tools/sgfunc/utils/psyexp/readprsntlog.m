function tbl = readprsntlog(fn)
% tbl = readprsntlog(fn)
%
% In Presentation logs, Unfortunately, # of delimiters across rows are
% INCONSISTENT! (kill me now)
%
% (cc) sgKIM, 2020.

% HARD CODED VARIABLES FOR PRESENTATION
HEADERLINES = 5;
VARIABLENAMES = {'Subject','Trial','EventType','Code','Time','TTime',...
  'TimeUncertainty','Duration','DurationUncertainty','ReqTime','ReqDur'};

%%
fid = fopen(fn,'r');
D = textscan(fid,'%s','HeaderLines',HEADERLINES ,'MultipleDelimsAsOne',1, ...
  'Delimiter',' ');
D = D{1};
fclose(fid);
nlines = size(D,1);
nvars = numel(VARIABLENAMES);
cellarray = cell(nlines,nvars);
X = table();
for iline = 1:nlines
  a = textscan(D{iline}, repmat('%q',[1 nvars]));
  x = [];
  for icol = 1:nvars
    if isempty(a{icol})
      x.(['x',num2str(icol)]) = nan;
    else
      x.(['x',num2str(icol)]) = string(a{icol}{1});
    end
  end
  x = struct2table(x);
  X = [X; x];
end
X.Properties.VariableNames = VARIABLENAMES;
fn_tmp = [tempname,'.csv'];
writetable(X,fn_tmp);
tbl = readtable(fn_tmp);
end