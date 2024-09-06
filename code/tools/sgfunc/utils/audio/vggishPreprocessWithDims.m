function [features, T_sec] = vggishPreprocessWithDims(audioIn,fs,varargin)
%VGGISHPREPROCESS Preprocess audio for VGGish feature extraction 
%   [features, TimeImage_sec] = vggishPreprocess(audioIn,fs) generates mel spectrograms from
%   the audio input, audioIn, that can be fed to the VGGish pretrained
%   network. fs is the sampling rate, in Hz. features is returned as a
%   96-by-64-by-1-by-N array, where 96 is the number of time hops, 64 is
%   the number of mel bands, and N is the number of generated mel
%   spectrograms.
%
%   features = vggishPreprocess(audioIn,fs,OverlapPercentage=OP) specifies
%   the overlap percentage between consecutive spectrograms. Specify OP as
%   a scalar in the range [0,100). If unspecified, OP defaults to 50.
%
%   Example:
%       % Use VGGish to extract embeddings from an audio signal.
%
%       % Read audio signal.
%       [audioIn,fs] = audioread("Counting-16-44p1-mono-15secs.wav");
%
%       % Preprocess the audio signal to extract mel spectrograms.
%       S = vggishPreprocess(audioIn,fs);
%
%       % Create VGGish network. This requires Deep Learning Toolbox.
%       net = vggish;
%
%       % Extract embeddings.
%       embeddings = predict(net,S);
%
% See also VGGISHEMBEDDINGS, VGGISH, YAMNET, CLASSIFYSOUND, YAMNETGRAPH

% Copyright 2020-2021 The MathWorks, Inc.

%#codegen

% Modified by dr.seunggoo.kim@gmail.com (2023-05-24)
% >>> NOT TO BE DISTRIBUTED PUBLICLY! <<<

narginchk(2, 6);

persistent filterBank numBands FFTLength FFTLengthOneSided FreqHz
if isempty(filterBank)
    FFTLength = 512;
    FFTLengthOneSided = 257;
    numBands = 64;
    frequencyRange = [125 7500];
    [filterBank, FreqHz] = designAuditoryFilterBank(...
      16e3, 'FrequencyScale','mel', 'FFTLength',512,...
      'FrequencyRange',frequencyRange, 'NumBands',64,...
      'Normalization','none', 'FilterBankDesignDomain','warped');
end

if isempty(coder.target)
    [fs, varargin{:}] = gather(fs, varargin{:});
end

fs0 = 16e3;

validateRequiredInputs(audioIn, fs);
params =  audio.internal.vggishPreprocessValidator(varargin{:});

if fs0 ~= fs
    x = resample(double(audioIn),fs0, double(fs));
    x = single(x);
else
    x = single(audioIn);
end

signalLength = size(x,1);
windowLength = 0.025 * fs0;
hopLength = 0.01 * fs0;
windowDuration = 0.96;
hopeDuration = 0.010;
featuresRate = 1.0 / hopeDuration;
frameLength = round(windowDuration * featuresRate);
minSignalLength = (frameLength-1)*hopLength + windowLength;
coder.internal.errorIf(signalLength<minSignalLength, ...
    'audio:vggish:SignalTooShort');

c = size(audioIn,2);
if isempty(coder.target) && isa(x,'gpuArray') 
    Y = audio.internal.buffer(x,windowLength,hopLength);
    Z   = abs(fft(Y.*hann(windowLength,'periodic'),FFTLength));
    binHigh = floor(FFTLength/2 + 1);
    Y       = head(Z,binHigh);
    b = size(Y,2)/c;
    T = audio.internal.buffer(1:numel(x), windowLength, hopLength);
else 
    [Y,~,T] = stft(x,'Window',hann(windowLength,'periodic'),...
        'OverlapLength',windowLength-hopLength,...
        'FFTLength',FFTLength,...
        'FrequencyRange','onesided');
    Y = abs(Y);
    [~,b,~] = size(Y);
    Y = reshape(Y,FFTLengthOneSided,b*c);
end

Z = filterBank * Y; % 257 filters -> 64 filters
Z = log(Z + single(params.LogAdditiveElement));
Z = permute(Z,[2 1]);

frameHopDuration = (100-params.OverlapPercentage) * windowDuration / 100;
frameHopLength = max(round(frameHopDuration * featuresRate), 1);

numHops = c * (floor((size(Z,1)/c-frameLength)/frameHopLength) + 1);
features = zeros(frameLength,numBands,1,numHops, 'like', x);
for index=1:c
    Z2 = audio.internal.buffer(Z((index-1)*b+1:index*b,:),frameLength,frameHopLength);
    L = size(Z2,2);
    M = L/numBands;
    Z2 = reshape(Z2,[size(Z2,1) M numBands]);
    Z2 = permute(Z2,[1 3 2]);
    Z2 = reshape(Z2,[size(Z2,1) numBands 1 M]);
    features(:,:, :,(index-1)*numHops/c + 1:index*numHops/c ) = Z2;
end

%F_Hz = FreqHz;
T2 = audio.internal.buffer(T,frameLength,frameHopLength);
T_sec = median(T2)/fs0;
assert(numel(T_sec) == size(features,4))
%assert(numel(F_Hz) == size(features,2))

fprintf('I: AudioIn length = %f sec\n',size(audioIn,1)/fs)
fprintf('I: Spectrogram dimensions = [%i time-bins, %i freq-bins, %i kerns, %i images]\n', ...
  size(features))
fprintf('I: Image center-times = [%.4f, %.4f] sec, dt = %.4f sec\n', ...
  T_sec([1 end]), diff(T_sec([1 2])) )

end

% -------------------------------------------------------------------------
% Validate required inputs
% -------------------------------------------------------------------------
function validateRequiredInputs(x,fs)
validateattributes(x,{'single','double'},...
    {'nonempty','2d','finite','real'}, ...
    'vggishPreprocess','audioIn')

validateattributes(fs,{'single','double'}, ...
    {'nonempty','positive','real','scalar','nonnan','finite'}, ...
    'vggishPreprocess','fs');
end
