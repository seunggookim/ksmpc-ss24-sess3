function LineLength = logthis(varargin)
%LOGTHIS FPRINTF with a header "[CallerName|TimeString-ISO8601] "
%
% LineLength = logthis('{FORMAT}', Inputs, ..., Verbosity={0|1|)
%
% (cc) 2021, sgKIM.


if numel(varargin)>1 && strcmpi(varargin{end-1}, 'verbosity') 
  if not(varargin{end})
    return
  else
    varargin(end-1:end) = [];
  end
end

DateStrNow = char(datetime('now'),'yyyy-MM-dd''T''HH:mm:ss');
% DateStrNow = char(datetime('now'),'yyyy-MM-dd HH:mm:ss');
st = dbstack;
callername = st(2).name;
LineLength = fprintf('[%s|%s] %s', callername, DateStrNow, sprintf(varargin{:}));
if ~nargout
  clear LineLength
end
end
