%DECODE_LMT Summary of this function goes here
%   Detailed explanation goes here
function decode_LMT_bootstrap
workingDir = pwd;

load([workingDir filesep 'ms.mat'])
load([workingDir filesep 'behav.mat'])

numBootstraps = 50;
Fs = 30;
min_speed_threshold = 5;
z_threshold = 2;

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
    if LMT_state(frame_i-1) == 4 && LMT_state(frame_i) == 1
        cum_distance(frame_i) = 0;
        cum_time(frame_i) = 0;
    else    
    cum_distance(frame_i,1) = cum_distance(frame_i-1)+abs(travelled_dist(frame_i));
    cum_time(frame_i,1) = cum_time(frame_i-1)+elapsed_time(frame_i);
    end 
end

single_trajectory_vec = interp_behav_vec;
single_trajectory_vec(LMT_state==2)=single_trajectory_vec(LMT_state==2)+134;
single_trajectory_vec(LMT_state==3)=single_trajectory_vec(LMT_state==3)+2*134;
single_trajectory_vec(LMT_state==4)=single_trajectory_vec(LMT_state==4)+3*134;

velocity = extract_velocity(interp_behav_vec, ca_time);
running_ts = velocity > min_speed_threshold;

%% Binarize
binarized_data = 0*ca_data;
for cell_i = 1:numNeurons
    binarized_data(:,cell_i) = extract_binary(ca_data(:,cell_i), Fs, z_threshold);
end

%% Defining training set
training_set_creation_method = 'random'; % 'odd', odd timestamps; 'first_portion', first portion of the recording; 3, 'random' random frames
training_set_portion = 0.5; % Portion of the recording used to train the decoder for method 2 and 3

training_ts = create_training_set(ca_time, training_set_creation_method, training_set_portion);
training_ts(running_ts == 0) = 0; % Exclude periods of immobility from the traing set

decoding_ts = ~training_ts; % Training timestamps are excluded
decoding_ts(running_ts == 0) = 0; % Periods of immobility are excluded

for bootstrap_i = 1:numBootstraps
    %% Sampling random epochs
    training_ts = create_training_set(ca_time, training_set_creation_method, training_set_portion);
    training_ts(running_ts == 0) = 0; % Exclude periods of immobility from the traing set
    decoding_ts = ~training_ts; % Training timestamps are excluded
    decoding_ts(running_ts == 0) = 0; % Periods of immobility are excludedcell_used = logical(zeros(size(ca_data,2),1));
    
    %% Sampling random cells
    num_cell_to_use = 100;
    cell_used = zeros(numNeurons,1);
    while sum(cell_used) < num_cell_to_use
        randomidx = ceil(rand*(length(cell_used)-1));
        cell_used(randomidx) = 1;
    end
    cell_used=logical(cell_used);
    
    %% Decode  
    [~, ~, mean_spatial_decoding_error(bootstrap_i), spatial_decoding_agreement(bootstrap_i)] = decode_1D(binarized_data(:,cell_used), interp_behav_vec, spatial_bin_vector, cell_used, training_ts, decoding_ts);
    [~, ~, mean_spatial_state_decoding_error(bootstrap_i), spatial_state_decoding_agreement(bootstrap_i)] = decode_1D(binarized_data(:,cell_used), single_trajectory_vec, cum_spatial_bin_vector, cell_used, training_ts, decoding_ts);
    [~, ~, mean_distance_decoding_error(bootstrap_i), distance_decoding_agreement(bootstrap_i)] = decode_1D(binarized_data(:,cell_used), cum_distance, cum_spatial_bin_vector, cell_used, training_ts, decoding_ts);
    [~, ~, mean_time_decoding_error(bootstrap_i), time_decoding_agreement(bootstrap_i)] = decode_1D(binarized_data(:,cell_used), cum_time, cum_temporal_bin_vector, cell_used, training_ts, decoding_ts);
    disp(['Done: ' num2str(bootstrap_i) '/' num2str(numBootstraps)])
end

%% Save the results
save([workingDir filesep 'LMT_decoding_results'],'mean_spatial_decoding_error','spatial_decoding_agreement','mean_spatial_state_decoding_error','spatial_state_decoding_agreement', 'mean_distance_decoding_error','distance_decoding_agreement','mean_time_decoding_error','time_decoding_agreement');

end

