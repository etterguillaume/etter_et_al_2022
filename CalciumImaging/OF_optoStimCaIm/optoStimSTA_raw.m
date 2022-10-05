function [STA_onset, STA_offset] = optoStimSTA_raw(ms, behav, cell_ID, OpticSigThreshold)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Parameters
Fs = round(1/mode(diff(ms.time/1000)));
dt = 1/Fs;
[bFilt,aFilt] = butter(2,  2/(30/2), 'low');

StimDuration = 5; % in s
StimDurationFrames = round(StimDuration*Fs);

%% Make sure time vectors contain unique time values
[unique_behavtime, IAbehav, ICbehav]=unique(behav.time);
[unique_mstime, IAms, ICms]=unique(ms.time);

calciumTrace = ms.RawTraces(IAms,cell_ID);

numFrames = length(calciumTrace);

%% Re-aligning calcium and behavior times
optoSignal = behav.optosignal;
InterpolatedOptoSignal = interp1(behav.time(IAbehav),optoSignal(IAbehav),ms.time(IAms));
smoothedOptoSignal = smooth(InterpolatedOptoSignal, 'sgolay');
normalizedOptoSignal = smoothedOptoSignal-min(smoothedOptoSignal);
normalizedOptoSignal = normalizedOptoSignal./max(normalizedOptoSignal);

%% Determine the stimulated epochs
%figure; plot(normalizedOptoSignal);
%[~,threshold] = ginput(1);

stimulatedEpochs = zeros(size(normalizedOptoSignal));
stimulatedEpochs(normalizedOptoSignal>OpticSigThreshold) = 1;

%% Determine stim onsets/offsets
stimDiff=diff(stimulatedEpochs);
stimDiff(end+1)=0;

stimOnset = stimDiff == 1;
stimOffset = stimDiff == -1;

stimOnset_idx = find(stimOnset == 1);
stimOffset_idx = find(stimOffset == 1);

%% STA analysis
%% Stim onset
STA_onset = zeros(length(stimOnset_idx),StimDurationFrames*2+1);

lines2remove = [];
for stim_i = 1:length(stimOnset_idx)
    if stimOnset_idx(stim_i) - StimDurationFrames <= 0 || stimOnset_idx(stim_i) + StimDurationFrames > length(stimOnset)
        lines2remove = stim_i;
    else      
        STA_onset(stim_i,:) = calciumTrace(stimOnset_idx(stim_i)-StimDurationFrames:stimOnset_idx(stim_i)+StimDurationFrames);
    end
end

% Remove incomplete baselines/stims
if ~isempty(lines2remove)
    STA_onset(lines2remove,:) = [];
end

%% Stim offset
STA_offset = zeros(length(stimOffset_idx),StimDurationFrames*2+1);

lines2remove = [];
for stim_i = 1:length(stimOffset_idx)
    if stimOffset_idx(stim_i) - StimDurationFrames <= 0 || stimOffset_idx(stim_i) + StimDurationFrames > length(stimOffset)
        lines2remove = stim_i;
    else      
        STA_offset(stim_i,:) = calciumTrace(stimOffset_idx(stim_i)-StimDurationFrames:stimOffset_idx(stim_i)+StimDurationFrames);
    end
end

% Remove incomplete baselines/stims
if ~isempty(lines2remove)
    STA_offset(lines2remove,:) = [];
end

end