function [zMI, p_values, baselineActivity, stimActivity, portBL] = optoStimMI(ms, behav, stimulatedEpochs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Parameters
Fs = round(1/mode(diff(ms.time/1000)));

%% Make sure time vectors contain unique time values
[unique_behavtime, IAbehav, ICbehav]=unique(behav.time);
[unique_mstime, IAms, ICms]=unique(ms.time);

numFrames = length(unique_mstime);

for cell_i = 1:size(ms.RawTraces,2)
    binarizedTrace = extract_binary(ms.RawTraces(IAms,cell_i), Fs,2);
    
%% Compute portion baseline activity
baselineActivity(cell_i) = sum(binarizedTrace == 1 & stimulatedEpochs == 0)/sum(stimulatedEpochs == 0);
stimActivity(cell_i) = sum(binarizedTrace == 1 & stimulatedEpochs == 1)/sum(stimulatedEpochs == 0);
portBL(cell_i) = stimActivity(cell_i)./baselineActivity(cell_i);

%% Compute probabilities and associated MI
MI = 0;
for activity_state = 0:1
    for stim_state = 0:1
    probActivityState = sum(binarizedTrace == activity_state)/numFrames;
    probStimState = sum(stimulatedEpochs == stim_state)/numFrames;
    jointProb = sum(binarizedTrace == activity_state & stimulatedEpochs == stim_state)/numFrames;

        if jointProb ~= 0
            MI = MI + jointProb*log2(jointProb./(probActivityState*probStimState));
        end
    end
end

%% Now let's shuffle and see if we obtain lower values
numShuffles = 1000;
shuffled_MI = zeros(numShuffles,1);
for shuffle_i = 1:numShuffles
    random_ts = ceil(rand*numFrames);
    shuffledBinarized = zeros(numFrames,1);
    % Permute the trace
    shuffledBinarized(1:random_ts) = binarizedTrace(end-random_ts+1:end);
    shuffledBinarized(random_ts+1:end) = binarizedTrace(1:end-random_ts);
    
    for activity_state = 0:1
        for stim_state = 0:1
        probActivityState = sum(shuffledBinarized == activity_state)/numFrames;
        probStimState = sum(stimulatedEpochs == stim_state)/numFrames;
        jointProb = sum(shuffledBinarized == activity_state & stimulatedEpochs == stim_state)/numFrames;

            if jointProb ~= 0
                shuffled_MI(shuffle_i) = shuffled_MI(shuffle_i) + jointProb*log2(jointProb./(probActivityState*probStimState));
            end
        end
    end
end

p_values(cell_i) = sum(shuffled_MI>MI)/numShuffles;
zMI(cell_i) = (MI - mean(shuffled_MI))/std(shuffled_MI);
end

end