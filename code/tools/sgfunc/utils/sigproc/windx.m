function x = windx(srate_hz,wind_ms,x,verbose)
% x = windx(srate_hz, wind_ms, x, [verbose])
%
% winding the begining and end of x for wind_ms
%
% INTPUS
%   srate_hz [1x1] sampling rate in Hz
%   wind_ms  [1x1] winding duration in ms
%   x        [NxP] signal in N samples in P channels (given N>P)
%   verbose  [1x1] 0=nothing, 1=text (default), or 2=text+figure
%
% OUTPUT
%   x        [NxP] wound signal in N samples in P channels
%
% (cc) 2019, sgKIM. solleo@gmail.com

%% Checking inputs
if ~exist('verbose','var')
  verbose = true;
end
d = size(x);
if isvector(x) % Introduced before R2006a
  x = reshape(x,[],1); % force it as a column-vector
else % if x is a matrix
  if d(1) < d(2) % if it's a matrix of row-vectors
    x = x';
  end
end

%% Lines from Tim's WIND()
[npts,nchns] = size(x);
wds = round(2*wind_ms/1000 * srate_hz);
if mod(wds,2)~=0
  wds = wds + 1;
end
t = linspace(-1*(pi/2), 1.5*pi, wds);
w = repmat((sin(t)'+1)/2,[1 nchns]);
idx = round(wds/2);
if verbose == 2 % for checking outputs
  x0 = x;
end
x(1:idx,:) = x(1:idx,:).*w(1:idx,:);
x(end-idx+1:end,:) = x(end-idx+1:end,:).*w(idx+1:wds,:);
if d(1) < d(2) % if it's a matrix of row-vectors
  x = x'; % putting back to the original dimension
end


%% Checking outputs
if verbose
  fprintf('#samples=%i, #channels=%i, samplingrate=%i Hz, wind=%i ms\n', ...
    npts, nchns, srate_hz, wind_ms);
end
if verbose == 2
  t0 = 0:1/srate_hz:(npts-1)/srate_hz;
  t1 = 1000*[0:1/srate_hz:(length(w)-1)/srate_hz];
  
  figure('position',[1 1 700 150]);
  subplot(231)
  plot(t0, x0)
  xlim([0 0.1]); ylabel('x'); xlabel('Time [s]');
  subplot(234)
  plot(t0, x0)
  xlim([t0(end)-0.1 t0(end)]);  ylabel('x'); xlabel('Time [s]')
  
  subplot(232)
  plot(t1(1:idx), w(1:idx,1))
  ylabel('w'); xlabel('Time [ms]')
  subplot(235)
  plot(t1(1:idx), w(end-idx+1:end,1))
  ylabel('w'); xlabel('Time [ms]')
  
  subplot(233)
  plot(t0, x)
  xlim([0 0.1]); ylabel('Wound x'); xlabel('Time [s]')
  subplot(236)
  plot(t0, x)
  xlim([t0(end)-0.1 t0(end)]);  ylabel('Wound x'); xlabel('Time [s]')
end

end
