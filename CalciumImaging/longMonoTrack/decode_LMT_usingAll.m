%DECODE_LMT Summary of this function goes here
%   Detailed explanation goes here
function decode_LMT_bootstrap
workingDir = pwd;

load([workingDir filesep 'ms.mat'])
load([workingDir filesep 'behav.mat'])
load([workingDir filesep 'LMT_results.mat'])

num_cell_to_use = 160;

numBootstraps = 50;
Fs = 30;
min_speed_threshold = 5;
z_threshold = 2;
lap_detect = 10;
decodeTempFilt = 2; % Temporal filter

temporalBinSize = 1;
spatialBinSize = 3;
max_cum_spatial_length = 536;
max_cum_temporal_length = 50;

spatial_bin_vector = 0:spatialBinSize:134+spatialBinSize;
spatial_bin_centers_vector = spatial_bin_vector + spatialBinSize/2;
spatial_bin_centers_vector(end) = [];

cum_spatial_bin_vector = 0:spatialBinSize:max_cum_spatial_length+spatialBinSize;
cum_spatial_bin_centers_vector = cum_spatial_bin_vector + spatialBinSize/2;
cum_spatial_bin_centers_vector(end) = [];

cum_temporal_bin_vector = 0:temporalBinSize:max_cum_temporal_length+temporalBinSize;
cum_temporal_bin_centers_vector = cum_temporal_bin_vector + temporalBinSize/2;
cum_temporal_bin_centers_vector(end) = [];

ca_time = ms.time/1000;
ca_data = ms.RawTraces;

behav_vec = behav.position(:,1);
behav_time = behav.time/1000;

optosignal = behav.optosignal;

%% Extract tone state
LMT_state = extractLMTState(optosignal);

[behav_time, IAbehav, ~] = unique(behav_time);
[ca_time, IAms, ~] = unique(ca_time);
ca_data = ca_data(IAms,:);
behav_vec = behav_vec(IAbehav,:);
LMT_state = LMT_state(IAbehav);
numNeurons = size(ca_data,2);

%% Interpolate behavior
interp_behav_vec = interpolate_behavior(behav_vec, behav_time, ca_time); % in the X dimension
interp_behav_vec(end) = interp_behav_vec(end-1);

%% Extract state
LMT_state = interp1(behav_time, LMT_state,ca_time,'nearest');

%% Extract traveled distance
travelled_dist = diff(interp_behav_vec);
travelled_dist(end+1)=0;

elapsed_time = diff(ca_time);
elapsed_time(end+1) = 0;

cum_distance(1) = 0;
cum_time(1) = 0;

for frame_i = 2:length(LMT_state)
    if LMT_state(frame_i) == 1 && interp_behav_vec(frame_i) < lap_detect
        cum_distance(frame_i) = 0;
        cum_time(frame_i) = 0;
    else
        cum_distance(frame_i) = cum_distance(frame_i-1)+abs(travelled_dist(frame_i));
        cum_time(frame_i) = cum_time(frame_i-1)+elapsed_time(frame_i);
    end
end

velocity = extract_velocity(interp_behav_vec, ca_time);
running_ts = velocity > min_speed_threshold;

%% Binarize
binarized_data = 0*ca_data;
for cell_i = 1:numNeurons
    binarized_data(:,cell_i) = extract_binary(ca_data(:,cell_i), Fs, z_threshold);
end

%% Defining training set
training_set_creation_method = 'random'; % 'odd', odd timestamps; 'first_portion', first portion of the recording; 3, 'random' random frames
training_set_portion = 0.9; % Portion of the recording used to train the decoder for method 2 and 3

for bootstrap_i = 1:numBootstraps
    %% Sampling random epochs
    training_ts = create_training_set(ca_time, training_set_creation_method, training_set_portion);
    training_ts(running_ts == 0) = 0; % Exclude periods of immobility from the traing set
    decoding_ts = ~training_ts; % Training timestamps are excluded
    decoding_ts(running_ts == 0) = 0; % Periods of immobility are excludedcell_used = logical(zeros(size(ca_data,2),1));
        
    %% Sampling random cells
    cell_used = zeros(numNeurons,1);
    while sum(cell_used) < num_cell_to_use
        randomidx = ceil(rand*(length(cell_used)-1));
        cell_used(randomidx) = 1;
    end
    cell_used=logical(cell_used);
    
    %% Create tuning curves
    for cell_i = 1:size(binarized_data,2)
        [~, ~, spatial_occupancy_vector, spatial_prob_being_active(cell_i), spatial_tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), interp_behav_vec, spatial_bin_vector, training_ts);
        [~, ~, distance_occupancy_vector, distance_prob_being_active(cell_i), distance_tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), cum_distance, cum_spatial_bin_vector, training_ts);
        [~, ~, time_occupancy_vector, time_prob_being_active(cell_i), time_tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), cum_time, cum_temporal_bin_vector, training_ts);
        
        spatial_tuning_curve_data(:,cell_i) = smoothPlaceFieldParam(spatial_tuning_curve_data(:,cell_i),5);
        time_tuning_curve_data(:,cell_i) = smoothPlaceFieldParam(time_tuning_curve_data(:,cell_i),10);
        distance_tuning_curve_data(:,cell_i) = smoothPlaceFieldParam(distance_tuning_curve_data(:,cell_i),15);
        
        random_ts = ceil(rand*length(ca_time));
        shuffled_binarized = circshift(binarized_data(:,cell_i),random_ts);
        
        [~, ~, ~, ~, shuffled_spatial_tuning_curve_data(:,cell_i)] = extract_1D_information(shuffled_binarized, interp_behav_vec, spatial_bin_vector, training_ts);
        [~, ~, ~, ~, shuffled_distance_tuning_curve_data(:,cell_i)] = extract_1D_information(shuffled_binarized, cum_distance, cum_spatial_bin_vector, training_ts);
        [~, ~, ~, ~, shuffled_time_tuning_curve_data(:,cell_i)] = extract_1D_information(shuffled_binarized, cum_time, cum_temporal_bin_vector, training_ts);
        
        shuffled_spatial_tuning_curve_data(:,cell_i) = smoothPlaceFieldParam(shuffled_spatial_tuning_curve_data(:,cell_i),5);
        shuffled_time_tuning_curve_data(:,cell_i) = smoothPlaceFieldParam(shuffled_time_tuning_curve_data(:,cell_i),10);
        shuffled_distance_tuning_curve_data(:,cell_i) = smoothPlaceFieldParam(shuffled_distance_tuning_curve_data(:,cell_i),15);
    end
    
    % Minimal a priori (use to remove experimental a priori)
    spatial_occupancy_vector = spatial_occupancy_vector./spatial_occupancy_vector*(1/length(spatial_occupancy_vector));
    distance_occupancy_vector = distance_occupancy_vector./distance_occupancy_vector*(1/length(distance_occupancy_vector));
    time_occupancy_vector = time_occupancy_vector./time_occupancy_vector*(1/length(time_occupancy_vector));
    
    %% Decode
    [spatial_decoded_probabilities] = bayesian_decode1D(binarized_data, spatial_occupancy_vector, spatial_prob_being_active, spatial_tuning_curve_data, cell_used);
    [spatial_decoded_probabilities] = bayesian_temporal_filter1D(spatial_decoded_probabilities,ca_time,decodeTempFilt);
    
    [distance_decoded_probabilities] = bayesian_decode1D(binarized_data, distance_occupancy_vector, distance_prob_being_active, distance_tuning_curve_data, cell_used);
    [distance_decoded_probabilities] = bayesian_temporal_filter1D(distance_decoded_probabilities,ca_time,decodeTempFilt);
    
    [time_decoded_probabilities] = bayesian_decode1D(binarized_data, time_occupancy_vector, time_prob_being_active, time_tuning_curve_data, cell_used);
    [time_decoded_probabilities] = bayesian_temporal_filter1D(time_decoded_probabilities,ca_time,decodeTempFilt);
    
    [shuffled_spatial_decoded_probabilities] = bayesian_decode1D(binarized_data, spatial_occupancy_vector, spatial_prob_being_active, shuffled_spatial_tuning_curve_data, cell_used);
    [shuffled_spatial_decoded_probabilities] = bayesian_temporal_filter1D(shuffled_spatial_decoded_probabilities,ca_time,decodeTempFilt);
    
    [shuffled_distance_decoded_probabilities] = bayesian_decode1D(binarized_data, distance_occupancy_vector, distance_prob_being_active, shuffled_distance_tuning_curve_data, cell_used);
    [shuffled_distance_decoded_probabilities] = bayesian_temporal_filter1D(shuffled_distance_decoded_probabilities,ca_time,decodeTempFilt);
    
    [shuffled_time_decoded_probabilities] = bayesian_decode1D(binarized_data, time_occupancy_vector, time_prob_being_active, shuffled_time_tuning_curve_data, cell_used);
    [shuffled_time_decoded_probabilities] = bayesian_temporal_filter1D(shuffled_time_decoded_probabilities,ca_time,decodeTempFilt);
    
    %% Let us now estimate the mouse location using the maximum a posteriori (MAP) value
    [spatial_max_decoded_prob, spatial_decoded_bin] = max(spatial_decoded_probabilities,[],1);
    decoded_position = spatial_bin_centers_vector(spatial_decoded_bin);
    [distance_max_decoded_prob, distance_decoded_bin] = max(distance_decoded_probabilities,[],1);
    decoded_distance = cum_spatial_bin_centers_vector(distance_decoded_bin);
    [time_max_decoded_prob, time_decoded_bin] = max(time_decoded_probabilities,[],1);
    decoded_time = cum_temporal_bin_centers_vector(time_decoded_bin);
    
    % same for shuffled surrogates
    [shuffled_spatial_max_decoded_prob, shuffled_spatial_decoded_bin] = max(shuffled_spatial_decoded_probabilities,[],1);
    shuffled_decoded_position = spatial_bin_centers_vector(shuffled_spatial_decoded_bin);
    [shuffled_distance_max_decoded_prob, shuffled_distance_decoded_bin] = max(shuffled_distance_decoded_probabilities,[],1);
    shuffled_decoded_distance = cum_spatial_bin_centers_vector(shuffled_distance_decoded_bin);
    [shuffled_time_max_decoded_prob, shuffled_time_decoded_bin] = max(shuffled_time_decoded_probabilities,[],1);
    shuffled_decoded_time = cum_temporal_bin_centers_vector(shuffled_time_decoded_bin);
    
    %% Discretize output variables
    % Location
    spatial_actual_bin = nan*interp_behav_vec;
    actual_position = nan*interp_behav_vec;
    for bin_i = 1:length(spatial_bin_vector)-1
        position_idx = find(interp_behav_vec>spatial_bin_vector(bin_i) & interp_behav_vec < spatial_bin_vector(bin_i+1));
        spatial_actual_bin(position_idx) = bin_i;
        actual_position(position_idx) = spatial_bin_centers_vector(bin_i);
    end
    
    % Distance
    distance_actual_bin = nan*cum_distance;
    actual_distance = nan*cum_distance;
    for bin_i = 1:length(cum_spatial_bin_vector)-1
        position_idx = find(cum_distance>cum_spatial_bin_vector(bin_i) & cum_distance < cum_spatial_bin_vector(bin_i+1));
        distance_actual_bin(position_idx) = bin_i;
        actual_distance(position_idx) = cum_spatial_bin_centers_vector(bin_i);
    end
    
    % Time
    time_actual_bin = nan*cum_time;
    actual_time = nan*cum_time;
    for bin_i = 1:length(cum_temporal_bin_vector)-1
        position_idx = find(cum_time>cum_temporal_bin_vector(bin_i) & cum_time < cum_temporal_bin_vector(bin_i+1));
        time_actual_bin(position_idx) = bin_i;
        actual_time(position_idx) = cum_temporal_bin_centers_vector(bin_i);
    end

    %% Compute errors
    spatial_decoded_bin(~decoding_ts) = nan;
    distance_decoded_bin(~decoding_ts) = nan;
    time_decoded_bin(~decoding_ts) = nan;
    shuffled_spatial_decoded_bin(~decoding_ts) = nan;
    shuffled_distance_decoded_bin(~decoding_ts) = nan;
    shuffled_time_decoded_bin(~decoding_ts) = nan;
    
    decoded_position(~decoding_ts) = nan;
    decoded_distance(~decoding_ts) = nan;
    decoded_time(~decoding_ts) = nan;
    shuffled_decoded_position(~decoding_ts) = nan;
    shuffled_decoded_distance(~decoding_ts) = nan;
    shuffled_decoded_time(~decoding_ts) = nan;

    spatial_actual_bin(~decoding_ts) = nan;
    actual_position(~decoding_ts) = nan;
    spatial_actual_bin = spatial_actual_bin';
    actual_position =  actual_position';
    
    distance_actual_bin(~decoding_ts) = nan;
    actual_distance(~decoding_ts) = nan;
    distance_actual_bin = distance_actual_bin';
    %actual_distance =  actual_distance';
    
    time_actual_bin(~decoding_ts) = nan;
    actual_time(~decoding_ts) = nan;
    time_actual_bin = time_actual_bin';
    %actual_time =  actual_time';
    
    %% Compute decoding agreement
%    decoding_agreement_vector = double(decoded_bin == actual_bin);
%    decoding_agreement_vector(isnan(decoded_bin)) = nan;
%    decoding_agreement_vector(isnan(actual_bin)) = nan;
%    decoding_agreement_vector(isnan(decoding_agreement_vector)) = [];
%    decoding_agreement = sum(decoding_agreement_vector)./length(decoding_agreement_vector);
    
    %% Compute decoding error
    spatial_decoding_error(bootstrap_i) = mean(abs(actual_position - decoded_position), 'omitnan');
    distance_decoding_error(bootstrap_i) = mean(abs(actual_distance - decoded_distance), 'omitnan');
    time_decoding_error(bootstrap_i) = mean(abs(actual_time - decoded_time), 'omitnan');
    
    shuffled_spatial_decoding_error(bootstrap_i) = mean(abs(actual_position - shuffled_decoded_position), 'omitnan');
    shuffled_distance_decoding_error(bootstrap_i) = mean(abs(actual_distance - shuffled_decoded_distance), 'omitnan');
    shuffled_time_decoding_error(bootstrap_i) = mean(abs(actual_time - shuffled_decoded_time), 'omitnan');
    
    disp(['Done: ' num2str(bootstrap_i) '/' num2str(numBootstraps)])
end

%% Save the results
save([workingDir filesep 'LMT_decoding_using_All_2gauss_160cells_2sFilt'],'spatial_decoding_error','distance_decoding_error','time_decoding_error','shuffled_spatial_decoding_error', 'shuffled_distance_decoding_error','shuffled_time_decoding_error');

end

