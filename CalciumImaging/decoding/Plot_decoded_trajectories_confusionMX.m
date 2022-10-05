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
decodeTempFilt = 2;

%% Interpolate
[interp_behav_vec] = interpolate_behavior(behav_vec, behav_time, ca_time);
interp_behav_vec(end) = interp_behav_vec(end-1);
[velocity] = extract_velocity(interp_behav_vec, ca_time);
min_speed_threshold = 5; % 2 cm.s-1
running_ts = velocity > min_speed_threshold;

%% Bin here
bin_size = 3;
bin_vector = min(interp_behav_vec):bin_size:max(interp_behav_vec)+bin_size; % start : bin_size : end
bin_centers_vector = bin_vector + bin_size/2;
bin_centers_vector(end) = [];

binarized_data = zeros(size(ca_data));
for cell_i = 1:size(ca_data,2)
    binarized_data(:,cell_i) = extract_binary(ca_data(:,cell_i), sampling_frequency, z_threshold);
end

training_set_creation_method = 'random'; % 'odd', odd timestamps; 'first_portion', first portion of the recording; 3, 'random' random frames
training_set_portion = 0.9; % Portion of the recording used to train the decoder for method 2 and 3

training_ts = create_training_set(ca_time, training_set_creation_method, training_set_portion);
training_ts(running_ts == 0) = 0; % Exclude periods of immobility from the traing set

decoding_ts = ~training_ts; % Training timestamps are excluded
decoding_ts(running_ts == 0) = 0; % Periods of immobility are excluded

for cell_i = 1:size(binarized_data,2)
    [~, ~, occupancy_vector, prob_being_active(cell_i), tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), interp_behav_vec, bin_vector, training_ts);
    tuning_curve_data(:,cell_i) = smoothPlaceFieldParam(tuning_curve_data(:,cell_i),5);
end

%% Decoding
occupancy_vector = occupancy_vector./occupancy_vector*(1/length(occupancy_vector));
cell_used = logical(ones(size(ca_data,2),1)); % Let us use every cell for now
[decoded_probabilities] = bayesian_decode1D(binarized_data, occupancy_vector, prob_being_active, tuning_curve_data, cell_used);
[decoded_probabilities] = bayesian_temporal_filter1D(decoded_probabilities,ca_time,decodeTempFilt);

%% Maximum a posteriori (MAP)
[max_decoded_prob, decoded_bin] = max(decoded_probabilities,[],1);
decoded_position = bin_centers_vector(decoded_bin);

%% Actual position
actual_bin = nan*interp_behav_vec;
actual_position = nan*interp_behav_vec;
for bin_i = 1:length(bin_vector)-1
    position_idx = find(interp_behav_vec>bin_vector(bin_i) & interp_behav_vec < bin_vector(bin_i+1));
    actual_bin(position_idx) = bin_i;
    actual_position(position_idx) = bin_centers_vector(bin_i);
end

%% Remove epochs of immobility
decoded_bin(~decoding_ts) = nan;
decoded_position(~decoding_ts) = nan;
decoded_probabilities(:,~decoding_ts) = nan;

actual_bin(~decoding_ts) = nan;
actual_position(~decoding_ts) = nan;
actual_bin = actual_bin';
actual_position =  actual_position';

%% Compute confusion matrix
confusion_matrix = zeros(length(bin_centers_vector),length(bin_centers_vector));

for actual_i = 1:length(bin_centers_vector)
   for decoded_i = 1:length(bin_centers_vector)
       confusion_matrix(actual_i,decoded_i) = sum(decoded_bin == decoded_i & actual_bin == actual_i)./sum(actual_bin == actual_i);
   end
end

% Plot the confusion matrix
figure
imagesc(bin_centers_vector, bin_centers_vector, confusion_matrix)
colormap viridis
title 'Confusion matrix'
xlabel 'Actual position (cm)'
ylabel 'Decoded position (cm)'