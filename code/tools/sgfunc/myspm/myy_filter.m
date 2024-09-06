function y = myy_filter(y, TR, Band)
% AllVolume = myy_filter(y, TR, Band)
%
% y: 2-D data [time x voxels] to filter
% TR: (sec)should be in seconds
% Band: [lowfreq,highfreq] (Hz)
%
% runs y_IdealFilter
% when the given volume is greater than 91*109*91*100, segments input into
% 20 segments for memory 
%
% (cc) 2015, sgKIM.  solleo@gmail.com   https://ggooo.wordpress.com

demean=0;
y = denan(y);
if numel(y) > 91*109*91*50
CUTNUMBER = 20;
% Filtering segmented for saving memory
fprintf('\n#Filtering...');
SegmentLength = ceil(size(y,2) / CUTNUMBER);
for iCut=1:CUTNUMBER
if iCut~=CUTNUMBER
Segment = (iCut-1)*SegmentLength+1 : iCut*SegmentLength;
else
Segment = (iCut-1)*SegmentLength+1 : size(y,2);
end
y(:,Segment) = y_IdealFilter(y(:,Segment), TR, Band, demean);
fprintf('.');
end
fprintf('\n');
else
y = y_IdealFilter(y, TR, Band, demean);
end
end

function [Data_Filtered] = y_IdealFilter(Data, SamplePeriod, Band, demean)
% FORMAT    [Data_Filtered] = y_IdealFilter(Data, SamplePeriod, Band, demean)
% Input:
% 	Data		    -	2D data matrix (nDimTimePoints * nTimeSeries)
% 	SamplePeriod	-   Sample period, i.e., 1/sample frequency. E.g., TR
%   Band            -   The frequency for filtering, 1*2 Array. Could be:
%                   [LowCutoff_HighPass HighCutoff_LowPass]: band pass filtering
%                   [0 HighCutoff_LowPass]: low pass filtering
%                   [LowCutoff_HighPass 0]: high pass filtering
% Output:
%	Data_Filtered       -   The data after filtering
%-----------------------------------------------------------
% Written by YAN Chao-Gan 120504 based on REST.
% The Nathan Kline Institute for Psychiatric Research, 140 Old Orangeburg Road, Orangeburg, NY 10962, USA
% Child Mind Institute, 445 Park Avenue, New York, NY 10022, USA
% The Phyllis Green and Randolph Cowen Institute for Pediatric Neuroscience, New York University Child Study Center, New York, NY 10016, USA
% ycg.yan@gmail.com


sampleFreq 	 = 1/SamplePeriod;
sampleLength = size(Data,1);
paddedLength = 2^nextpow2(sampleLength);
LowCutoff_HighPass = Band(1);
HighCutoff_LowPass = Band(2);

% Get the frequency index
if (LowCutoff_HighPass >= sampleFreq/2) % All high stop
idxLowCutoff_HighPass = paddedLength/2 + 1;
else % high pass, such as freq > 0.01 Hz
idxLowCutoff_HighPass = ceil(LowCutoff_HighPass * paddedLength * SamplePeriod + 1);
end

if (HighCutoff_LowPass>=sampleFreq/2)||(HighCutoff_LowPass==0) % All low pass
idxHighCutoff_LowPass = paddedLength/2 + 1;
else % Low pass, such as freq < 0.08 Hz
idxHighCutoff_LowPass = fix(HighCutoff_LowPass * paddedLength * SamplePeriod + 1);
end

FrequencyMask = zeros(paddedLength,1);
FrequencyMask(idxLowCutoff_HighPass:idxHighCutoff_LowPass,1) = 1;
FrequencyMask(paddedLength-idxLowCutoff_HighPass+2:-1:paddedLength-idxHighCutoff_LowPass+2,1) = 1;

%FrequencySetZero_Index = find(FrequencyMask==0);

%Remove the mean before zero padding for boundary artifacts
DataMean = repmat(mean(Data),size(Data,1),1);
Data = Data - DataMean;
Data = [Data;zeros(paddedLength -sampleLength,size(Data,2))]; %padded with zero
Data = fft(Data);
Data(FrequencyMask==0,:) = 0; % sharp transition would cause ringing artifact?!
Data = ifft(Data);
Data_Filtered = Data(1:sampleLength,:);
if exist('demean','var')&&(demean==0)
Data_Filtered = Data_Filtered + DataMean;
end
end
