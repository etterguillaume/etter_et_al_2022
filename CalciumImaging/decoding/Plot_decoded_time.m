load('behav.mat')
load('ms.mat')

ca_data = ms.RawTraces;
ca_time = ms.time/1000;
behav_time=behav.time/1000;
behav_vec = behav.position(:,1);

[behav_time, IAbehav, ~] = unique(behav_time);
[ca_time, IAms, ~] = unique(ca_time);
ca_data = ca_data(IAms,:);
behav_vec = behav_vec(IAbehav,:);
numNeurons = size(ca_data,2);
sampling_frequency = 30; % This data set has been sampled at 30 images per second
z_threshold = 2;

%% Interpolate
[interp_behav_vec] = interpolate_behavior(behav_vec, behav_time, ca_time);
interp_behav_vec(end) = interp_behav_vec(end-1);
[velocity] = extract_velocity(interp_behav_vec, ca_time);
min_speed_threshold = 5; % 2 cm.s-1
running_ts = velocity > min_speed_threshold;

%% Bin time here
temporalBinSize = 1;
max_cum_temporal_length = 50;
cum_temporal_bin_vector = 0:temporalBinSize:max_cum_temporal_length+temporalBinSize;
cum_temporal_bin_centers_vector = cum_temporal_bin_vector + temporalBinSize/2;
cum_temporal_bin_centers_vector(end) = [];

elapsed_time = diff(ca_time);
elapsed_time(end+1) = 0;

%% Extract opto state
optosignal = behav.optosignal;
LMT_state = extractLMTState(optosignal);
LMT_state = LMT_state(IAbehav);
LMT_state = interp1(behav_time, LMT_state,ca_time,'nearest');

lap_detect = 10;
for frame_i = 2:length(LMT_state)
    if LMT_state(frame_i) == 1 && interp_behav_vec(frame_i) < lap_detect
        cum_time(frame_i) = 0;
    else
        cum_time(frame_i) = cum_time(frame_i-1)+elapsed_time(frame_i);
    end
end

binarized_data = zeros(size(ca_data));
for cell_i = 1:size(ca_data,2)
    binarized_data(:,cell_i) = extract_binary(ca_data(:,cell_i), sampling_frequency, z_threshold);
end

%% Extract tuning curves
for cell_i = 1:size(binarized_data,2)
    [~, ~, occupancy_vector, prob_being_active(cell_i), tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), cum_time, cum_temporal_bin_vector, running_ts);
end

%% Decoding
occupancy_vector = occupancy_vector./occupancy_vector*(1/length(occupancy_vector));
cell_used = logical(ones(size(ca_data,2),1)); % Let us use every cell for now
[decoded_probabilities] = bayesian_decode1D(binarized_data, occupancy_vector, prob_being_active, tuning_curve_data, cell_used);

%% Maximum a posteriori (MAP)
[max_decoded_prob, decoded_bin] = max(decoded_probabilities,[],1);
decoded_time = cum_temporal_bin_centers_vector(decoded_bin);

%% Actual time
actual_bin = nan*cum_time;
actual_time = nan*cum_time;
for bin_i = 1:length(cum_temporal_bin_vector)-1
    time_idx = find(cum_time>cum_temporal_bin_vector(bin_i) & cum_time < cum_temporal_bin_vector(bin_i+1));
    actual_bin(time_idx) = bin_i;
    actual_time(time_idx) = cum_temporal_bin_centers_vector(bin_i);
end

%% Remove epochs of immobility
decoded_bin(~running_ts) = nan;
decoded_time(~running_ts) = nan;
%decoded_probabilities(:,~running_ts) = nan;

%% Plot actual vs decoded time
figure
subplot(2,1,1)
imagesc(ca_time,cum_temporal_bin_centers_vector,decoded_probabilities)
colormap Viridis
title 'Posterior probabilities'
xlabel 'Time (s)'
ylabel 'Elapsed time (s)'
ax1 = gca;
ax1.CLim = [0 0.1];
ax1.XLim = [175 230];
ax1.YDir = 'normal';
subplot(2,1,2)
plot(ca_time,actual_time)
hold on
plot(ca_time, decoded_time)
title 'Actual versus decoded position'
xlabel 'Time (s)'
ylabel 'Location on the track (cm)'
ax2 = gca;
ax2.XLim = [175 230]; % Let's plot a single trajectory
linkaxes([ax1 ax2], 'x')