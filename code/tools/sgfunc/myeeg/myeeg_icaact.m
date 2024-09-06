function EEG = myeeg_icaact(EEG)
% EEG = myeeg_icaact(EEG)
%
% NOTE: EEG.icaweight is a demixing matrix W such that
% A = XSW
%
% where X is a measured signal [#times x #chans], 
%       S is a whitening matrix [#chans x #comp1], 
%       W is a demixing matrix [#comp1 x #comp2], 
%   and A is a soucre activity [#times x #comp2].
% #comp1 equals to #chans when X is full rank otherwise # of eigenvariates
%   with suprathreshold eigenvalues.
% #comp2 is given by users.
%
% (cc) sgKIM, 2019.
EEG.icaact = eeg_getdatact(EEG, 'component', 1:size(EEG.icaweights,1));
end