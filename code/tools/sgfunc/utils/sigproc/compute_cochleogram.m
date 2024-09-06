function [ccg,times,freqs] = compute_cochleogram(audio,fs,frameSize_sec)
% COCHLEOGRAM computes a cochleogram (gammatone spectrogram) using MATLAB 
% Audio Toolbox. This script is simply a wrapper of an example script
% of GAMMATONEFILTERBANK
%
% USAGES
% [ccg, t_sec, f_Hz] = cochleogram(audiosignal, sampleRate_Hz)
%
% [ccg, t_sec, f_Hz] = cochleogram(audiosignal, gammatoneFilterBankObject)
% 
% [ccg, t_sec, f_Hz] = cochleogram(..., frameSize_sec)
% 
%
% (cc?) 2019, sgKIM.
%
% SEE ALSO: GAMMATONEFILTERBANK

if strcmp(class(fs), 'gammatoneFilterBank')
  gammaFiltBank = fs;
  fs = get(gammaFiltBank,'SampleRate');
else
  gammaFiltBank = gammatoneFilterBank('SampleRate',fs);
end
audioOut = gammaFiltBank(audio);

% Computes the energy-per-band using 50 ms windows + 25 ms overlap
if ~exist('frameSize_sec','var')
  frameSize_sec = 0.050;
end
samplesPerFrame = round(frameSize_sec*fs);
samplesOverlap = round(frameSize_sec*0.5*fs);
buff = dsp.AsyncBuffer(numel(audio));
write(buff, audioOut.^2);
sink = dsp.SignalSink;
while buff.NumUnreadSamples > 0
    currentFrame = read(buff, samplesPerFrame, samplesOverlap);
    sink(mean(currentFrame,1))
end

% Convert the energy values to dB
ccg = 20*log10(sink.Buffer'+eps);
times = ((samplesPerFrame-samplesOverlap)/fs)*(0:size(ccg,2)-1);
freqs = getCenterFrequencies(gammaFiltBank);

if 0 % for TEST
  figure
  mypcolor(times,freqs,ccg)
  xlabel('Time [s]'); ylabel('Freq [Hz]')
  cb = colorbar;
  ylabel(cb,'Energy [dB/Hz]')
end
end