function [randX, rnd_theta] = randomize_phase(X, rnd_theta, keep_cov)
% [randX, rnd_theta] = randomize_phase(X, rnd_theta, keep_cov)
%
% X must be a column vector or a matrix of column vectors

% REF: https://www.mathworks.com/matlabcentral/answers/451578-fft-and-ifft-random-phases#answer_379416

% Check for NaN
assert(all(~isnan(X(:))), '[%s] %i/%i Not-a-Numbers!: fix it as you please and come back!', ...
  mfilename, sum(isnan(X(:))), numel(X))

% Set indices for odd/even-sample signals
[N,P] = size(X);
if mod(N,2) % odd-sample signal
  nAngles = (N-1)/2;
  IdxPos = 2:(N+1)/2;
  IdxNeg = (N+1)/2+1:N;
else % even-sample signal
  nAngles = N/2-1;
  IdxPos = 2:N/2;
  IdxNeg = N/2+2:N;
end

% Take out the mean
X0 = mean(X,1);
X = X - X0;

% 1. Get spectrum
Y = fft(X); 

% 2. Add random phase shifts (negative for conjugates), preserve DC offset
if ~exist('rnd_theta','var') || isempty(rnd_theta)
  if ~exist('keep_cov','var')
    keep_cov = false;
  end
  if keep_cov
    rnd_theta = -pi + repmat((2*pi).*rand(nAngles,1), [1 P]);
  else
    rnd_theta = -pi + (2*pi).*rand(nAngles,P);
  end
end
Y(IdxPos,:) = Y(IdxPos,:).*exp(1i*rnd_theta);
Y(IdxNeg,:) = Y(IdxNeg,:).*exp(-1i*flipud(rnd_theta));

% 3. Reconstruct phase-randomized data
randX = ifft(Y);

% Put back the mean
randX = randX + X0;

end

function TEST()
%%
t = [0:0.0001:1]';
X = cos(2*pi*10*t) + cos(2*pi*30*t);
size(X)
Y = randomize_phase(X);
clf; 
plot(t, X, t, Y)
assert( rms(abs(fft(X))-abs(fft(Y))) < 1e-12, 'Quite BIG error???')
%%
end

