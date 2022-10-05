function extract_stimSpeedActivityScrambled
%EXTRACT_STIMSPEEDACTIVITY Summary of this function goes here
%   Detailed explanation goes here

%% Load data
workingDir = pwd;

load([workingDir filesep 'ms.mat']);
load([workingDir filesep 'behav.mat']);

%% Parameters
sampling_frequency = 30;
z_threshold = 2;
binSize = 2;%cm.s-1

speed_bin_vector = [0 5 30];
speed_bins_centers = speed_bin_vector+binSize/2;

OpticSigThreshold = 0.25;
Fs = round(1/mode(diff(ms.time/1000)));
[bFilt,aFilt] = butter(2,  2/(30/2), 'low');

optoSignal = behav.optosignal;
smoothedOptoSignal = smooth(optoSignal, 'sgolay');
normalizedOptoSignal = smoothedOptoSignal-min(smoothedOptoSignal);
normalizedOptoSignal = normalizedOptoSignal./max(normalizedOptoSignal);

figure; plot(normalizedOptoSignal);

%% Make sure time vectors contain unique time values
[behav_time, IAbehav, ICbehav]=unique(behav.time/1000);
[ca_time, IAms, ICms]=unique(ms.time/1000);

optoSignal = behav.optosignal;

%% Interpolate signals
InterpolatedOptoSignal = interp1(behav_time,optoSignal(IAbehav),ca_time);
[interp_behav_vec] = interpolate_behavior(behav.position(IAbehav,:), behav_time, ca_time);
interp_behav_vec(end) = interp_behav_vec(end-1);

%% Detect stimulation epochs
smoothedOptoSignal = smooth(InterpolatedOptoSignal, 'sgolay');
normalizedOptoSignal = smoothedOptoSignal-min(smoothedOptoSignal);
normalizedOptoSignal = normalizedOptoSignal./max(normalizedOptoSignal);

stimulatedEpochs = zeros(size(normalizedOptoSignal));
stimulatedEpochs(normalizedOptoSignal>OpticSigThreshold) = 1;
stimulatedEpochs=logical(stimulatedEpochs);

%% Extract velocity
[velocity] = extract_velocity(interp_behav_vec, ca_time);

%% Binarize data and compute relationship to speed
numNeurons = size(ms.RawTraces,2);
for cell_i = 1:numNeurons
    binarizedTrace = extract_binary(ms.RawTraces(:,cell_i), sampling_frequency, z_threshold);
    [BL_MI{cell_i}, ~, BL_occupancy, BL_activity{cell_i}, BL_speed_activity{cell_i}(:)] = extract_1D_information(binarizedTrace, velocity, speed_bin_vector, ~stimulatedEpochs);
    [stim_MI{cell_i}, ~, stim_occupancy, stim_activity{cell_i}, stim_speed_activity{cell_i}(:)] = extract_1D_information(binarizedTrace, velocity, speed_bin_vector, stimulatedEpochs);
end

%% Output
%save something here
save([workingDir filesep 'scrambledStim_speedActivity.mat'],'BL_MI','stim_MI','BL_occupancy','stim_occupancy','BL_activity','stim_activity','BL_speed_activity','stim_speed_activity');

end

