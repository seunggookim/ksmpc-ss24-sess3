function varargout = yamnetPreprocessWithDims(audioIn,fs,varargin)
% yamnetPreprocess Preprocess audio for YAMNet classification 
%   features = yamnetPreprocess(audioIn,fs) generates mel spectrograms from
%   the audio input, audioIn, that can be fed to the YAMNet pretrained
%   network. fs is the sampling rate, in Hz. features is returned as a 
%   96-by-64-by-1-by-N array, where N is the number of generated mel
%   spectrograms.
%
%   features = yamnetPreprocess(audioIn,fs,OverlapPercentage=OP) specifies
%   the overlap percentage between consecutive spectrograms. Specify OP as
%   a scalar in the range [0,100). If unspecified, OP defaults to 50.
%
%   Example:
%       % Use YAMNet to classify the sound in an audio signal.
%
%       % Read audio signal.
%       [audioIn,fs] = audioread("Counting-16-44p1-mono-15secs.wav");
%
%       % Preprocess the audio signal to extract mel spectrograms.
%       S = yamnetPreprocess(audioIn,fs);
%
%       % Create YAMNet network.
%       net = yamnet;
%
%       % Classify the spectrogram images.
%       classes = classify(net,S);
% 
%       % Classify the audio signal as the most frequently occurring sound.
%       mySound = mode(classes)
%
% See also YAMNET, VGGISHEMBEDDINGS, VGGISH, VGGISHPREPROCESS, 
%          CLASSIFYSOUND, YAMNETGRAPH

% Copyright 2020-2021 The MathWorks, Inc.

%#codegen

narginchk(2, 4);

[varargout{1:2}] = vggishPreprocessWithDims(audioIn,fs,'LogAdditiveElement',1e-3, varargin{:});
