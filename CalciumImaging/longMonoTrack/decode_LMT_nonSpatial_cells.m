%DECODE_LMT Summary of this function goes here
%   Detailed explanation goes here
function decode_LMT_spatial
workingDir = pwd;

load([workingDir filesep 'ms.mat'])
load([workingDir filesep 'behav.mat'])
load([workingDir filesep 'LMT_results.mat'])

num_cell_to_use = 32;

numBootstraps = 50;
Fs = 30;
min_speed_threshold = 5;
z_threshold = 2;
decodeTempFilt = 2; % Temporal filter
spatialBinSize = 3;

%% Load and identify place cells
cells2use = spatial_MI_pvalue > 0.05;

spatial_bin_vector = 0:spatialBinSize:134+spatialBinSize;
spatial_bin_centers_vector = spatial_bin_vector + spatialBinSize/2;
spatial_bin_centers_vector(end) = [];

ca_time = ms.time/1000;
ca_data = ms.RawTraces;

behav_vec = behav.position(:,1);
behav_time = behav.time/1000;

%% Extract tone state
[behav_time, IAbehav, ~] = unique(behav_time);
[ca_time, IAms, ~] = unique(ca_time);
ca_data = ca_data(IAms,:);
behav_vec = behav_vec(IAbehav,:);
numNeurons = size(ca_data,2);

%% Interpolate behavior
interp_behav_vec = interpolate_behavior(behav_vec, behav_time, ca_time); % in the X dimension
interp_behav_vec(end) = interp_behav_vec(end-1);

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
        num_cell_to_use = sum(cells2use);
        warning(['Only ' num2str(sum(cells2use)) ' significant place cells. Readjusting boostrap sample size accordingly...']);
    end

    cell_used = zeros(numNeurons,1);
    while sum(cell_used) < num_cell_to_use
        CellsIDX = find(cells2use==1);
        randomNum = ceil(rand*(length(CellsIDX)));
        randomidx= CellsIDX(randomNum);
        if cells2use(randomidx) == 1
            cell_used(randomidx) = 1;
        end
    end
    cell_used=logical(cell_used);
    
    %% Create tuning curves
    for cell_i = 1:size(binarized_data,2)
        [~, ~, spatial_occupancy_vector, spatial_prob_being_active(cell_i), spatial_tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), interp_behav_vec, spatial_bin_vector, training_ts);
        spatial_tuning_curve_data(:,cell_i) = smoothPlaceFieldParam(spatial_tuning_curve_data(:,cell_i),5);
        
        random_ts = ceil(rand*length(ca_time));
        shuffled_binarized = circshift(binarized_data(:,cell_i),random_ts);
        
        [~, ~, ~, ~, shuffled_spatial_tuning_curve_data(:,cell_i)] = extract_1D_information(shuffled_binarized, interp_behav_vec, spatial_bin_vector, training_ts);
        shuffled_spatial_tuning_curve_data(:,cell_i) = smoothPlaceFieldParam(shuffled_spatial_tuning_curve_data(:,cell_i),5);
    end
    
    % Minimal a priori (use to remove experimental a priori)
    spatial_occupancy_vector = spatial_occupancy_vector./spatial_occupancy_vector*(1/length(spatial_occupancy_vector));
    
    %% Decode
    [spatial_decoded_probabilities] = bayesian_decode1D(binarized_data, spatial_occupancy_vector, spatial_prob_being_active, spatial_tuning_curve_data, cell_used);
    [spatial_decoded_probabilities] = bayesian_temporal_filter1D(spatial_decoded_probabilities,ca_time,decodeTempFilt);
    
    [shuffled_spatial_decoded_probabilities] = bayesian_decode1D(binarized_data, spatial_occupancy_vector, spatial_prob_being_active, shuffled_spatial_tuning_curve_data, cell_used);
    [shuffled_spatial_decoded_probabilities] = bayesian_temporal_filter1D(shuffled_spatial_decoded_probabilities,ca_time,decodeTempFilt);
    
    %% Let us now estimate the mouse location using the maximum a posteriori (MAP) value
    [spatial_max_decoded_prob, spatial_decoded_bin] = max(spatial_decoded_probabilities,[],1);
    decoded_position = spatial_bin_centers_vector(spatial_decoded_bin);
    
    % same for shuffled surrogates
    [shuffled_spatial_max_decoded_prob, shuffled_spatial_decoded_bin] = max(shuffled_spatial_decoded_probabilities,[],1);
    shuffled_decoded_position = spatial_bin_centers_vector(shuffled_spatial_decoded_bin);
    
    %% Discretize output variables
    % Location
    spatial_actual_bin = nan*interp_behav_vec;
    actual_position = nan*interp_behav_vec;
    for bin_i = 1:length(spatial_bin_vector)-1
        position_idx = find(interp_behav_vec>spatial_bin_vector(bin_i) & interp_behav_vec < spatial_bin_vector(bin_i+1));
        spatial_actual_bin(position_idx) = bin_i;
        actual_position(position_idx) = spatial_bin_centers_vector(bin_i);
    end

    %% Compute errors
    spatial_decoded_bin(~decoding_ts) = nan;
    shuffled_spatial_decoded_bin(~decoding_ts) = nan;
    
    decoded_position(~decoding_ts) = nan;
    shuffled_decoded_position(~decoding_ts) = nan;

    spatial_actual_bin(~decoding_ts) = nan;
    actual_position(~decoding_ts) = nan;
    spatial_actual_bin = spatial_actual_bin';
    actual_position =  actual_position';
    
    %% Compute decoding error
    spatial_decoding_error(bootstrap_i) = mean(abs(actual_position - decoded_position), 'omitnan');    
    shuffled_spatial_decoding_error(bootstrap_i) = mean(abs(actual_position - shuffled_decoded_position), 'omitnan');

    disp(['Done: ' num2str(bootstrap_i) '/' num2str(numBootstraps)])
end

%% Save the results
save([workingDir filesep 'LMT_spatial_decoding_using_32nonSpatial_cells'],'spatial_decoding_error', 'shuffled_spatial_decoding_error');

end

