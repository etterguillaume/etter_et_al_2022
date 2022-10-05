%DECODE_LMT Summary of this function goes here
%   Detailed explanation goes here
function decode_LMT_distance_cells
workingDir = pwd;

load([workingDir filesep 'ms.mat'])
load([workingDir filesep 'behav.mat'])
load([workingDir filesep 'LMT_results.mat'])

num_cell_to_use = 25;

numBootstraps = 50;
Fs = 30;
min_speed_threshold = 5;
z_threshold = 2;
lap_detect = 10;
decodeTempFilt = 2; % Temporal filter

spatialBinSize = 3;
max_cum_spatial_length = 536;

%% Load and identify place cells
cells2use = distance_MI_pvalue < 0.05;

cum_spatial_bin_vector = 0:spatialBinSize:max_cum_spatial_length+spatialBinSize;
cum_spatial_bin_centers_vector = cum_spatial_bin_vector + spatialBinSize/2;
cum_spatial_bin_centers_vector(end) = [];

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

cum_distance(1) = 0;

for frame_i = 2:length(LMT_state)
    if LMT_state(frame_i) == 1 && interp_behav_vec(frame_i) < lap_detect
        cum_distance(frame_i) = 0;
    else
        cum_distance(frame_i) = cum_distance(frame_i-1)+abs(travelled_dist(frame_i));
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
    if sum(cells2use)<num_cell_to_use
        num_cell_to_use = sum(cells2use)
        warning(['Only ' num2str(sum(cells2use)) ' significant distance cells. Readjusting boostrap sample size accordingly...']);
    end
        
    cell_used = zeros(numNeurons,1);
    while sum(cell_used) < num_cell_to_use
        cellsIDX = find(cells2use==1);
        randomNum = ceil(rand*(length(cellsIDX)));
        randomidx= cellsIDX(randomNum);
        if cells2use(randomidx) == 1
            cell_used(randomidx) = 1;
        end
    end
    cell_used=logical(cell_used);
    
    %% Create tuning curves
    for cell_i = 1:size(binarized_data,2)
        [~, ~, distance_occupancy_vector, distance_prob_being_active(cell_i), distance_tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), cum_distance, cum_spatial_bin_vector, training_ts);
        
        distance_tuning_curve_data(:,cell_i) = smoothPlaceFieldParam(distance_tuning_curve_data(:,cell_i),15);
        
        random_ts = ceil(rand*length(ca_time));
        shuffled_binarized = circshift(binarized_data(:,cell_i),random_ts);
        
        [~, ~, ~, ~, shuffled_distance_tuning_curve_data(:,cell_i)] = extract_1D_information(shuffled_binarized, cum_distance, cum_spatial_bin_vector, training_ts);
        
        shuffled_distance_tuning_curve_data(:,cell_i) = smoothPlaceFieldParam(shuffled_distance_tuning_curve_data(:,cell_i),15);
        
    end
    
    % Minimal a priori (use to remove experimental a priori)
    distance_occupancy_vector = distance_occupancy_vector./distance_occupancy_vector*(1/length(distance_occupancy_vector));
    
    %% Decode
    [distance_decoded_probabilities] = bayesian_decode1D(binarized_data, distance_occupancy_vector, distance_prob_being_active, distance_tuning_curve_data, cell_used);
    [distance_decoded_probabilities] = bayesian_temporal_filter1D(distance_decoded_probabilities,ca_time,decodeTempFilt);
    
    [shuffled_distance_decoded_probabilities] = bayesian_decode1D(binarized_data, distance_occupancy_vector, distance_prob_being_active, shuffled_distance_tuning_curve_data, cell_used);
    [shuffled_distance_decoded_probabilities] = bayesian_temporal_filter1D(shuffled_distance_decoded_probabilities,ca_time,decodeTempFilt);
    
    %% Let us now estimate the mouse location using the maximum a posteriori (MAP) value
    [distance_max_decoded_prob, distance_decoded_bin] = max(distance_decoded_probabilities,[],1);
    decoded_distance = cum_spatial_bin_centers_vector(distance_decoded_bin);
    
    % same for shuffled surrogates
    [shuffled_distance_max_decoded_prob, shuffled_distance_decoded_bin] = max(shuffled_distance_decoded_probabilities,[],1);
    shuffled_decoded_distance = cum_spatial_bin_centers_vector(shuffled_distance_decoded_bin);
    
    %% Discretize output variables
    % Distance
    distance_actual_bin = nan*cum_distance;
    actual_distance = nan*cum_distance;
    for bin_i = 1:length(cum_spatial_bin_vector)-1
        position_idx = find(cum_distance>cum_spatial_bin_vector(bin_i) & cum_distance < cum_spatial_bin_vector(bin_i+1));
        distance_actual_bin(position_idx) = bin_i;
        actual_distance(position_idx) = cum_spatial_bin_centers_vector(bin_i);
    end

    %% Compute errors
    distance_decoded_bin(~decoding_ts) = nan;
    shuffled_distance_decoded_bin(~decoding_ts) = nan;
    
    decoded_distance(~decoding_ts) = nan;
    shuffled_decoded_distance(~decoding_ts) = nan;
    
    distance_actual_bin(~decoding_ts) = nan;
    actual_distance(~decoding_ts) = nan;
    distance_actual_bin = distance_actual_bin';
    %actual_distance =  actual_distance';
    
    %% Compute decoding error
    distance_decoding_error(bootstrap_i) = mean(abs(actual_distance - decoded_distance), 'omitnan');
    
    shuffled_distance_decoding_error(bootstrap_i) = mean(abs(actual_distance - shuffled_decoded_distance), 'omitnan');
    
    disp(['Done: ' num2str(bootstrap_i) '/' num2str(numBootstraps)])
end

%% Save the results
save([workingDir filesep 'LMT_decoding_using_25distance_cells'],'distance_decoding_error', 'shuffled_distance_decoding_error');

end

