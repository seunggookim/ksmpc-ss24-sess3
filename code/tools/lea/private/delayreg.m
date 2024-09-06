function X = delayreg(x, TimeMask, Delays, isZeropad)
% X = delayreg(x, TimeMask, Delays, isZeropad)
%
% INPUTS:
%   x          [numerical: #times x #features]
%   TimeMask   [logcial: #samples x 1]
%   Delays     [numerical: 1 x #delays>]
%   isZeropad  [logical: false (default) | true]
%
% OUTPUTS:
% For X = [x1, x2, x3, ...]
%     DELAYS = [d0, d1, d2,...],
% Returns  X = [x1(d0), x2(d0), x3(d0), ..., x1(d1), x2(d2), ...]
%
% because this is easier to concatenate along the columns for various
% feature sets: [X1 X2]
%
% without zeropadding:
%  - if x(delay) is NaNs or undefined, return an error 
%
% (cc) 2021-2022, sgKIM

% y(t) = [x(t-tau1), x(t-tau2), x(t-tau3), ...]*beta
%
% e.g., tau = [0, 1, 2, ...]    y(t) is modelled by previous x's
%       tau = [-2, -1, 0]       y(t) is modelled by following x's

if ~exist('isZeropad','var'), isZeropad = false; end


TimeMask = not(not(TimeMask));
[~,nFeats] = size(x);
nSamples = sum(TimeMask);
nLags = numel(Delays);
X = zeros(sum(TimeMask), nFeats*numel(Delays), class(x));
for i = 1:nLags
  Idx = delay(TimeMask, Delays(i));  % delaying a logical index array
  nZeros = nSamples - sum(Idx);      % count #zeros to pad
  
  if ~isZeropad
    assert(~nZeros, ['For given TimeMask & Delays, there is no valid ',...
      'sample of X. Try different settings or zero-padding.'])
  end
  
  X(:,i+(0:(nFeats-1))*nLags) = [...
    zeros(nZeros*(Delays(i)>0),nFeats)     % zeropadding for positive lags
    x(Idx,:)                               % slicing input
    zeros(nZeros*(Delays(i)<0),nFeats) ... % zeropadding for negative lags
    ];
end

% denan:
X(isnan(X)) = 0;

end

function U = delay(I, tau)
if (tau<=0) % a negative delay
  U = [false(-tau,1); I(1:end+tau)];
else % a positive delay
  U = [I(1+tau:end); false(tau,1)];
end
end
