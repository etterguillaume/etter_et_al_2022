%DECODE_LMT Summary of this function goes here
%   Detailed explanation goes here
function decode_LMT_temporal_cells
workingDir = pwd;

load([workingDir filesep 'ms.mat'])
load([workingDir filesep 'behav.mat'])
load([workingDir filesep 'LMT_results.mat'])

num_cell_to_use = 14;

numBootstraps = 50;
Fs = 30;
min_speed_threshold = 5;
z_threshold = 2;
lap_detect = 10;
decodeTempFilt = 2; % Temporal filter

temporalBinSize = 1;
spatialBinSize = 3;
max_cum_temporal_length = 50;

%% Load and identify place cells
cells2use = time_MI_pvalue < 0.05;

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
elapsed_time = diff(ca_time);
elapsed_time(end+1) = 0;

cum_time(1) = 0;

for frame_i = 2:length(LMT_state)
    if LMT_state(frame_i) == 1 && interp_behav_vec(frame_i) < lap_detect
        cum_time(frame_i) = 0;
    else
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
    if sum(cells2use)<num_cell_to_use
        num_cell_to_use = sum(cells2use)
        warning(['Only ' num2str(sum(cells2use)) ' significant time cells. Readjusting boostrap sample size accordingly...']);
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
        [~, ~, time_occupancy_vector, time_prob_being_active(cell_i), time_tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), cum_time, cum_temporal_bin_vector, training_ts);
        time_tuning_curve_data(:,cell_i) = smoothPlaceFieldParam(time_tuning_curve_data(:,cell_i),10);
        
        random_ts = ceil(rand*length(ca_time));
        shuffled_binarized = circshift(binarized_data(:,cell_i),random_ts);
        
        [~, ~, ~, ~, shuffled_time_tuning_curve_data(:,cell_i)] = extract_1D_information(shuffled_binarized, cum_time, cum_temporal_bin_vector, training_ts);
        shuffled_time_tuning_curve_data(:,cell_i) = smoothPlaceFieldParam(shuffled_time_tuning_curve_data(:,cell_i),10);
        
    end
    
    % Minimal a priori (use to remove experimental a priori)
    time_occupancy_vector = time_occupancy_vector./time_occupancy_vector*(1/length(time_occupancy_vector));
    
    %% Decode
    [time_decoded_probabilities] = bayesian_decode1D(binarized_data, time_occupancy_vector, time_prob_being_active, time_tuning_curve_data, cell_used);
    [time_decoded_probabilities] = bayesian_temporal_filter1D(time_decoded_probabilities,ca_time,decodeTempFilt);

    [shuffled_time_decoded_probabilities] = bayesian_decode1D(binarized_data, time_occupancy_vector, time_prob_being_active, shuffled_time_tuning_curve_data, cell_used);
    [shuffled_time_decoded_probabilities] = bayesian_temporal_filter1D(shuffled_time_decoded_probabilities,ca_time,decodeTempFilt);
    
    %% Let us now estimate the mouse location using the maximum a posteriori (MAP) value
    [time_max_decoded_prob, time_decoded_bin] = max(time_decoded_probabilities,[],1);
    decoded_time = cum_temporal_bin_centers_vector(time_decoded_bin);
    
    % same for shuffled surrogates
    [shuffled_time_max_decoded_prob, shuffled_time_decoded_bin] = max(shuffled_time_decoded_probabilities,[],1);
    shuffled_decoded_time = cum_temporal_bin_centers_vector(shuffled_time_decoded_bin);
    
    %% Discretize output variables
    % Time
    time_actual_bin = nan*cum_time;
    actual_time = nan*cum_time;
    for bin_i = 1:length(cum_temporal_bin_vector)-1
        position_idx = find(cum_time>cum_temporal_bin_vector(bin_i) & cum_time < cum_temporal_bin_vector(bin_i+1));
        time_actual_bin(position_idx) = bin_i;
        actual_time(position_idx) = cum_temporal_bin_centers_vector(bin_i);
    end

    %% Compute errors
    time_decoded_bin(~decoding_ts) = nan;
    shuffled_time_decoded_bin(~decoding_ts) = nan;
    
    decoded_time(~decoding_ts) = nan;
    shuffled_decoded_time(~decoding_ts) = nan;
    
    time_actual_bin(~decoding_ts) = nan;
    actual_time(~decoding_ts) = nan;
    time_actual_bin = time_actual_bin';
    %actual_time =  actual_time';

    %% Compute decoding error
    time_decoding_error(bootstrap_i) = mean(abs(actual_time - decoded_time), 'omitnan');
    
    shuffled_time_decoding_error(bootstrap_i) = mean(abs(actual_time - shuffled_decoded_time), 'omitnan');
    
    disp(['Done: ' num2str(bootstrap_i) '/' num2str(numBootstraps)])
end

%% Save the results
save([workingDir filesep 'LMT_decoding_using_14temporal_cells'],'time_decoding_error','shuffled_time_decoding_error');

end

