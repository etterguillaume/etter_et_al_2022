%% Assess coding of location, direction, distance, time, and state
%% Load data here
workingDir = pwd;

load([workingDir filesep 'ms.mat']);
load([workingDir filesep 'behav.mat']);

%% Parameters
numShuffles = 1000;
Fs = 30;
min_speed_threshold = 5;
z_threshold = 2;
lap_detect = 10;

temporalBinSize = 1;
spatialBinSize = 3;
max_cum_spatial_length = 536;
max_cum_temporal_length = 10;

spatial_bin_vector = 0:spatialBinSize:134+spatialBinSize;
spatial_bin_centers_vector = spatial_bin_vector + spatialBinSize/2;
spatial_bin_centers_vector(end) = [];

cum_spatial_bin_vector = 0:spatialBinSize:max_cum_spatial_length+spatialBinSize;
cum_spatial_bin_centers_vector = cum_spatial_bin_vector + spatialBinSize/2;
cum_spatial_bin_centers_vector(end) = [];

cum_temporal_bin_vector = 0:temporalBinSize:max_cum_temporal_length+temporalBinSize;
cum_temporal_bin_centers_vector = cum_temporal_bin_vector + temporalBinSize/2;
cum_temporal_bin_centers_vector(end) = [];

%% Extract variables from data structure
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

%single_trajectory_vec = interp_behav_vec;
%single_trajectory_vec(LMT_state==2)=single_trajectory_vec(LMT_state==2)+134;
%single_trajectory_vec(LMT_state==3)=single_trajectory_vec(LMT_state==3)+2*134;
%single_trajectory_vec(LMT_state==4)=single_trajectory_vec(LMT_state==4)+3*134;

%% Extract velocity
velocity = extract_velocity(interp_behav_vec, ca_time);
running_ts = velocity > min_speed_threshold;

%% Binarize
binarized_data = 0*ca_data;
for cell_i = 1:numNeurons
    binarized_data(:,cell_i) = extract_binary(ca_data(:,cell_i), Fs, z_threshold);
end

%% Extract tuning curves
for cell_i = 1:numNeurons
    [spatial_MI(cell_i), ~, spatial_occupancy_vector, spatial_prob_being_active(cell_i), spatial_tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), interp_behav_vec, spatial_bin_vector, running_ts);
    %[spatial_state_MI, ~, spatial_state_occupancy_vector, spatial_state_prob_being_active(cell_i), spatial_state_tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), single_trajectory_vec, cum_spatial_bin_vector, running_ts);
    [distance_MI(cell_i), ~, distance_occupancy_vector, distance_prob_being_active(cell_i), distance_tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), cum_distance, cum_spatial_bin_vector, running_ts);
    [time_MI(cell_i), ~, time_occupancy_vector, time_prob_being_active(cell_i), time_tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), cum_time, cum_temporal_bin_vector, running_ts);
    
    for shuffle_i = 1:numShuffles
    random_ts = ceil(rand*length(ca_time));
    shuffled_binarized = circshift(binarized_data(:,cell_i),random_ts); 
    
    [shuffled_spatial_MI(shuffle_i), ~, ~, ~, ~] = extract_1D_information(shuffled_binarized, interp_behav_vec, spatial_bin_vector, running_ts);
    %[shuffled_spatial_state_MI(shuffle_i), ~, ~, ~, ~] = extract_1D_information(shuffled_binarized, single_trajectory_vec, cum_spatial_bin_vector, running_ts);
    [shuffled_distance_MI(shuffle_i), ~, ~, ~, ~] = extract_1D_information(shuffled_binarized, cum_distance, cum_spatial_bin_vector, running_ts);
    [shuffled_time_MI(shuffle_i), ~, ~, ~, ~] = extract_1D_information(shuffled_binarized, cum_time, cum_temporal_bin_vector, running_ts);
    end
    
    spatial_MI_pvalue(cell_i) = sum(shuffled_spatial_MI>spatial_MI(cell_i))/numShuffles;
    %spatial_state_MI_pvalue(cell_i) = sum(shuffled_spatial_state_MI>spatial_state_MI)/numShuffles;
    distance_MI_pvalue(cell_i) = sum(shuffled_distance_MI>distance_MI(cell_i))/numShuffles;
    time_MI_pvalue(cell_i) = sum(shuffled_time_MI>time_MI(cell_i))/numShuffles;
    
    spatial_zMI(cell_i) = (spatial_MI(cell_i)-mean(shuffled_spatial_MI))/std(shuffled_spatial_MI);
    %spatial_state_zMI(cell_i) = (spatial_state_MI-mean(shuffled_spatial_state_MI))/std(shuffled_spatial_state_MI);
    distance_zMI(cell_i) = (distance_MI(cell_i)-mean(shuffled_distance_MI))/std(shuffled_distance_MI);
    time_zMI(cell_i) = (time_MI(cell_i)-mean(shuffled_time_MI))/std(shuffled_time_MI);
    
    disp([num2str(round(cell_i/numNeurons*100)) '% completed'])
end

save([workingDir filesep 'LMT_results_10s.mat'],'time_MI', 'distance_MI', 'spatial_MI', 'spatial_zMI', 'distance_zMI', 'time_zMI', 'time_MI_pvalue', 'distance_MI_pvalue', 'spatial_MI_pvalue', 'spatial_tuning_curve_data', 'distance_tuning_curve_data', 'time_tuning_curve_data');

