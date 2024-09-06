function itpc = myft_compute_itpc(F)
% computes inter-trial phase coherence (ITPC)
%
% itpc = myft_compute_itpc(F)
%
% F is a complex Fourier coefficient in #trials x #chan x #freq x #time
% (FieldTrip convention)
%
% REF1: http://www.fieldtriptoolbox.org/faq/itc/
% REF2: Delorme A, Makeig S. EEGLAB: an open source toolbox for analysis of single-trial EEG dynamics including independent component analysis. J Neurosci Methods. 2004 Mar 15;134(1):9-21. pdf
% REF3: Catherine Tallon-Baudry, Olivier Bertrand, Claude Delpuech and Jacques Pernier, Journal of Neuroscience 1 July 1996, 16 (13) 4240-4249; DOI: https://doi.org/10.1523/JNEUROSCI.16-13-04240.1996
%{
[REF3]
The phase-locking of the oscillatory burst can be evaluated in the
time–frequency domain by adapting the “phase-averaging” methods previously proposed in the frequency domain by Jervis et al. (1983). We
consider the normalized complex time-varying energy of each single trial i: 
Pi(t, f0) = w(t, f0).si(t)/|w(t,f0).si(t)|. Averaging these quantities across
single trials leads to a complex value related to the phase distribution of
each time–frequency region around t and f0. The modulus of this value
will be called the “phase-locking factor.” It ranges from 0 (purely nonphase-locked activity) to 1 (strictly phase-locked activity)
%}

itpc = squeeze( abs( sum(F./abs(F),1) ) / size(F,1) );
end