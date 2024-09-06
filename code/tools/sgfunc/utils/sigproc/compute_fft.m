function [F,freq_hz] = compute_fft(y,fs)
% [F,freq_hz] = compute_fft(y,fs)
% only returns nonnegative frequency components
% (cc) 2019, sgKIM.
F = fft(y);
num_tp = size(y,1);
freq_hz = [0:ceil(num_tp/2)-1]/(num_tp/fs); % positive frequency bins
F = F(1:numel(freq_hz),:);
end