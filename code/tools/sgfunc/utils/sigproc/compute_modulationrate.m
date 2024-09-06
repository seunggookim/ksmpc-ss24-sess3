function [power_dB,freq_hz] = compute_modulationrate(y,fs)
% (cc) 2019, sgKIM.

H = hilbert(y);
E = abs(H); % magnitude of the analytic signal is the envelope
F = fft(E);
power_dB = 20.*log10(abs(F)); % 20*log10(magnitude) = 10*log10(mag^2)
freq_hz = [0:length(y)/2-1]/(length(y)/fs); % positive frequency bins
power_dB = power_dB(1:numel(freq_hz),:); % discarding negative frequencies
power_dB = power_dB - max(power_dB(:));
% if exist('fig','var') && fig
%   t = [0:length(y)-1]/fs;
%   
%   figure
%   subplot(221)
%   H0 = hilbert(y0);
%   E0 = sqrt(real(H0).^2 + imag(H0).^2);
%   plot(t,y0,'b', t,E0,'r')
%   xlim([1 1.5])
%   title('Unfiltered')
%   
%   subplot(222)
%   plot(freq_hz, power_dB)
%   xlim([0 100])
%   
%   subplot(223)
%   plot(t,y,'b', t,E,'r')
%   xlim([1 1.5])
%   title(['Filtered at [',num2str(bpf),'] Hz'])
%   
%   subplot(224)
%   F0 = fft(E0);
%   plot(freq_hz, 20*log10(abs(F(1:numel(freq_hz))+eps)))
%   xlim([0 100])
% end
end