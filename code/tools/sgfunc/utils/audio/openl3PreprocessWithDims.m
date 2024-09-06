function [features, times_s] = openl3PreprocessWithDims(audioIn,fs,varargin)
%OPENL3PREPROCESS Preprocess audio for OpenL3 feature extraction
%   [features,freqs_hz,times_s] = openl3PreprocessWithDims(audioIn,fs) 
%   generates spectrograms from the audio input, audioIn, that can be fed 
%   to the OpenL3 pretrained network. fs is the sampling rate, in Hz.
%
%   [features,...] = openl3PreprocessWithDims(audioIn,fs,OverlapPercentage=OP) 
%   specifies
%   the overlap percentage between consecutive spectrograms. Specify OP as a
%   scalar in the range [0,100). The default is 90.
%
%   [features,...] = openl3PreprocessWithDims(audioIn,fs,SpectrumType=SType) 
%   specifies
%   the type of output spectrum type, as "mel128", "mel256", or "linear".
%   The default is "mel128".
%
%   Example:
%       % Use OpenL3 to extract embeddings from an audio signal.
%
%       % Read audio signal.
%       [audioIn,fs] = audioread("Counting-16-44p1-mono-15secs.wav");
%
%       % Preprocess the audio signal to extract mel spectrograms.
%       S = openl3Preprocess(audioIn,fs);
%
%       % Create OpenL3 network. This requires Deep Learning Toolbox.
%       net = openl3;
%
%       % Extract embeddings.
%       embeddings = predict(net,S);
%
% See also OPENL3, OPENL3EMBEDDINGS, VGGISHEMBEDDINGS, VGGISH

% Copyright 2020-2021 The MathWorks, Inc.

%#codegen

% Modified by dr.seunggoo.kim@gmail.com (2023-05-24)

narginchk(2, 6);

if isempty(coder.target)
    [fs, varargin{:}] = gather(fs, varargin{:});
end

validateRequiredInputs(audioIn, fs);
params =  audio.internal.openl3PreprocessValidator(varargin{:});
isLinear = params.IsLinear;
fs0 = 48e3;

persistent wRealMel wImagMel wRealLin wImagLin fb128 fb256

if isempty(fb128) && strcmp(params.SpectrumType,'mel128')
    [fb128] = designAuditoryFilterBank(fs0,'FFTLength', 2048, 'NumBands', 128,...
        'Normalization', 'bandwidth',...
        'OneSided', true);   
end
if isempty(fb256)&& strcmp(params.SpectrumType,'mel256')
    [fb256] = designAuditoryFilterBank(fs0,'FFTLength', 2048, 'NumBands', 256,...
        'Normalization', 'bandwidth',...
        'OneSided', true);
end

if ~isLinear && isempty(wRealMel)
    win = single(hann(2048, 'periodic'));
    N = 2048;
    n = (0:N-1)';
    wRealMel = zeros(1025,2048,'single');
    wImagMel = zeros(1025,2048,'single');
    for indx = 1:1025
        wRealMel(indx,:) = -sin(2*pi*(indx-1)*n/N).*win;
        wImagMel(indx,:) = cos(2*pi*(indx-1)*n/N).*win;
    end
end

if isLinear && isempty(wRealLin)  
    win = single(hann(512, 'periodic'));
    N = 512;
    n = (0:N-1).';
    wRealLin = zeros(257,512,'single');
    wImagLin = zeros(257,512,'single');
    for indx = 1:257
        wRealLin(indx,:) = -sin(2*pi*(indx-1)*n/N).*win;
        wImagLin(indx,:) = cos(2*pi*(indx-1)*n/N).*win;
    end
end

if fs0 ~= fs
    x = resample(double(audioIn), fs0, double(fs));
    x = single(x);
else
    x = single(audioIn);
end

frameLength = fs0;
hopLength  = round(fs0 - fs0 * params.OverlapPercentage/100);
frames = audio.internal.buffer(x, frameLength, hopLength);
T = audio.internal.buffer((1:numel(x))', frameLength, hopLength);

SpectrumType = params.SpectrumType;
if isLinear
    wReal = wRealLin;
    wImag = wImagLin;
    NFFT = 512;
    NFFTHalf = 257;
    sz = [NFFTHalf 197 1];
    pad = 0;
else
    wReal = wRealMel;
    wImag = wImagMel;
    NFFT = 2048;
    NFFTHalf = 1025;
    if strcmp(SpectrumType, 'mel128')
        melSize = 128;
    else
        melSize = 256;
    end
    sz = [melSize 199 1];
    pad = (242*(sz(2)-1) + NFFT - 48000)/2;
end

f = audio.internal.buffer([...
    zeros(pad, size(frames, 2)); ...
    frames; ...
    zeros(pad, size(frames,2))], NFFT, 242);
T = audio.internal.buffer(...
    [-inf(pad, size(T,2)); T; inf(pad, size(T,2))], NFFT, 242);

S0 = (wReal * f).^2 + (wImag * f).^2;
S = reshape(S0, [NFFTHalf sz(2) 1 size(S0,2)/sz(2)]);
T = reshape(median(T),[sz(2) size(S0,2)/sz(2)]);

if ~isLinear
    S = reshape(S, NFFTHalf, []);
    if strcmp(SpectrumType, 'mel128')
        fb = fb128;
    else
       fb = fb256;
    end
    Sout = fb * S;
    Sout = reshape(Sout, melSize, 199, 1, []);
else
    Sout = S;
end

Sout = sqrt(Sout);
Sout = 10 * log10(max(Sout , single(1e-10)));
M = max(Sout, [], [1 2 3]);
Sout = bsxfun(@minus, Sout, M);
features = max(Sout, -80);

times_s  = median(T)/fs0;

fprintf('I: AudioIn length = %f sec\n',size(audioIn,1)/fs)
fprintf('I: Spectrogram dimensions = [%i freq-bins, %i time-bins, %i kerns, %i images]\n', size(features))
fprintf('I: Image center-times = [%.4f, %.4f] sec, dt = %.4f sec\n', times_s([1 end]), diff(times_s([1 2])) )

end

% -------------------------------------------------------------------------
% Validate required inputs
% -------------------------------------------------------------------------
function validateRequiredInputs(x,fs)
validateattributes(x,{'single','double'},...
    {'nonempty','2d','finite','real'}, ...
    'openl3Preprocess','audioIn')

validateattributes(fs,{'single','double'}, ...
    {'nonempty','positive','real','scalar','nonnan','finite'}, ...
    'openl3Preprocess','fs');
end
