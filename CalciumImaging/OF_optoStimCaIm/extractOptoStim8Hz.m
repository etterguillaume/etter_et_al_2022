function [stimulatedEpochs, OpticSigThreshold] = extractOptoStim8Hz(ms,behav)
%OPTOSTIMEXTRACTSIGTHRESHOLD Summary of this function goes here
%   Detailed explanation goes here

%% Parameters
Fs = round(1/mode(diff(ms.time/1000)));
stimDuration = 5;
frameDuration = stimDuration*Fs;
frameDuration = round(0.8*frameDuration); % This prevents potential misdetection of 

%% Set threshold
optoSignal = behav.optosignal;
normalizedOptoSignal = optoSignal-min(optoSignal);
normalizedOptoSignal = normalizedOptoSignal./max(normalizedOptoSignal);

figure; plot(normalizedOptoSignal);
[~,OpticSigThreshold] = ginput(1);
%OpticSigThreshold = 0.2;

%% Make sure time vectors contain unique time values
[unique_behavtime, IAbehav, ICbehav]=unique(behav.time);
[unique_mstime, IAms, ICms]=unique(ms.time);

%% Re-aligning calcium and behavior times
InterpolatedOptoSignal = interp1(behav.time(IAbehav),optoSignal(IAbehav),ms.time(IAms));
normalizedOptoSignal = InterpolatedOptoSignal-min(InterpolatedOptoSignal);
normalizedOptoSignal = normalizedOptoSignal./max(normalizedOptoSignal);

stimulatedEpochs = zeros(size(normalizedOptoSignal));
stimulatedEpochs(normalizedOptoSignal>OpticSigThreshold) = 1;

tempStimulatedEpochs=0*stimulatedEpochs;

stimOn = 0;

for frame_i = 2:length(stimulatedEpochs)-frameDuration
    if stimOn == 0
        if stimulatedEpochs(frame_i-1)==0 & stimulatedEpochs(frame_i)==1 & sum(stimulatedEpochs(frame_i:frame_i+frameDuration-1)) >=4
            stimOnsets(frame_i) = 1;
            stimOn = 1;
        end
    else
        if stimulatedEpochs(frame_i-1)==1 & stimulatedEpochs(frame_i)==0 & sum(stimulatedEpochs(frame_i:frame_i+frameDuration-1)) == 0
            stimOffsets(frame_i) = 1;
            stimOn = 0;
        end
    end
    tempStimulatedEpochs(frame_i) = stimOn;
end

stimulatedEpochs = tempStimulatedEpochs;

end

