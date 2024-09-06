function s = create_ripple(t, f, tmf, smf, phi)
% creates a ripple function (i.e., 2D cosine function)
%
% USAGE
% s = create_ripple(t, f, tmf, smf, phi)
%
% INPUTS
% t   [1x#TimeBins] time vector defining the function in sec
% f   [1x#FreqBins] frequency vector defining the function in Hz
% tmf [1x1] temporal modulation frequency in Hz
% smf [1x1] spectral modulation frequency in cycle/Hz
% phi [1x1] initial phase at (0,0) in degrees
%
% OUTPUTS
% s   [#TimeBins x #FreqBins]
%
% REF: Chi et al., 1999
% (cc) 2019, sgKIM

s = zeros(numel(t), numel(f));
for ti=1:numel(t)
  for fi=1:numel(f)
    s(ti,fi) = cos(2*pi*tmf*t(ti) + 2*pi*smf*f(fi) + phi/180*pi);
  end
end
end
