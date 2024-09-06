function [ctpc,ctsc] = myft_compute_ctpc (F)
% computes cross-trial phase coherence (CTPC) and cross-trial squared
% amplitude coherenece (CTSC)
%
% [ctpc, ctsc] = compute_ctpc (F)
%
% F is a complex Fourier coefficient in #trials x #chans x #freqs x #times
% (FieldTrip convention)
% Ref: Lou & Poeppel, 2007, Neuron, doi:10.1016/j.neuron.2007.06.004.
%
% for each freq f and time t:
% CTPC = mean_over_trial(cos(theta))^2 + mean_over_trial(sin(theta))^2
% CTSC = sqrt(mean_over_trial((amp^2 - mean(amp^2))^2)) / mean(amp^2)

theta = angle(F);
ctpc = squeeze( nanmean(cos(theta),1).^2 + nanmean(sin(theta),1).^2 );
SA = abs(F).^2;
ctsc = squeeze( sqrt(nanmean((SA.^2 - nanmean(SA.^2)).^2)) ./ nanmean(SA.^2) );
end

