function [power_dB,freq_hz] = compute_power(y,fs)
% [power_dB,freq_hz] = compute_power(y,fs)
% (cc) 2019, sgKIM.
F = fft(y);
power_dB = 20*log10(abs(F)+eps); % 20*log10(magnitude) = 10*log10(mag^2)
freq_hz = [0:length(y)/2-1]/(length(y)/fs); % positive frequency bins
power_dB = power_dB(1:numel(freq_hz)); % discarding negative frequencies

end