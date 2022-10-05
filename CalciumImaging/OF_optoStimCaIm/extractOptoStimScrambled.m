function [stimulatedEpochs,OpticSigThreshold] = extractOptoStimScrambled(ms,behav)
%OPTOSTIMEXTRACTSIGTHRESHOLD Summary of this function goes here
%   Detailed explanation goes here

%% Parameters
Fs = round(1/mode(diff(ms.time/1000)));
[bFilt,aFilt] = butter(2,  2/(30/2), 'low');

%% Set threshold
optoSignal = behav.optosignal;
smoothedOptoSignal = smooth(optoSignal, 'sgolay');
normalizedOptoSignal = smoothedOptoSignal-min(smoothedOptoSignal);
normalizedOptoSignal = normalizedOptoSignal./max(normalizedOptoSignal);

figure; plot(normalizedOptoSignal);
[~,OpticSigThreshold] = ginput(1);

%% Make sure time vectors contain unique time values
[unique_behavtime, IAbehav, ICbehav]=unique(behav.time);
[unique_mstime, IAms, ICms]=unique(ms.time);

%% Re-aligning calcium and behavior times
optoSignal = behav.optosignal;
InterpolatedOptoSignal = interp1(behav.time(IAbehav),optoSignal(IAbehav),ms.time(IAms));
smoothedOptoSignal = smooth(InterpolatedOptoSignal, 'sgolay');
normalizedOptoSignal = smoothedOptoSignal-min(smoothedOptoSignal);
normalizedOptoSignal = normalizedOptoSignal./max(normalizedOptoSignal);

stimulatedEpochs = zeros(size(normalizedOptoSignal));
stimulatedEpochs(normalizedOptoSignal>OpticSigThreshold) = 1;

end

