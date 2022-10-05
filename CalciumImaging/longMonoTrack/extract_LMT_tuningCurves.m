%% Assess coding of location, direction, distance, time, and state
%% Load data here
workingDir = pwd;

load([workingDir filesep 'ms.mat']);
load([workingDir filesep 'behav.mat']);

%% Parameters
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
    if LMT_state(frame_i-1) == 4 && LMT_state(frame_i) == 1
        cum_distance(frame_i) = 0;
        cum_time(frame_i) = 0;
    else    
    cum_distance(frame_i) = cum_distance(frame_i-1)+abs(travelled_dist(frame_i));
    cum_time(frame_i) = cum_time(frame_i-1)+elapsed_time(frame_i);
    end 
end

%% Create a single trajectory vector that spans across 4 tone states
% z_delta_x = diff(interp_behav_vec);
% z_delta_x(end+1) = 0;
% z_delta_x(isnan(z_delta_x)) = 0;
% z_delta_x = zscore(z_delta_x);
% 
% right_trajectories = interp_behav_vec;
% right_trajectories(z_delta_x<0.2) = NaN;
% 
% left_trajectories = interp_behav_vec;
% left_trajectories(z_delta_x>-0.2) = NaN;
% 
% left_trajectories1 = 2*max(interp_behav_vec) - left_trajectories;
% right_trajectories2 = 2*max(interp_behav_vec) + right_trajectories;
% left_trajectories2 = 4*max(interp_behav_vec) - left_trajectories;

% single_trajectory_vec = right_trajectories*nan;
% single_trajectory_vec(LMT_state==1) = right_trajectories(LMT_state==1);
% single_trajectory_vec(LMT_state==2) = left_trajectories1(LMT_state==2);
% single_trajectory_vec(LMT_state==3) = right_trajectories2(LMT_state==3);
% single_trajectory_vec(LMT_state==4) = left_trajectories2(LMT_state==4);

single_trajectory_vec = interp_behav_vec;
single_trajectory_vec(LMT_state==2)=single_trajectory_vec(LMT_state==2)+134;
single_trajectory_vec(LMT_state==3)=single_trajectory_vec(LMT_state==3)+2*134;
single_trajectory_vec(LMT_state==4)=single_trajectory_vec(LMT_state==4)+3*134;

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
    [actual_spatial_MI(cell_i), ~, spatial_occupancy_vector, spatial_prob_being_active(cell_i), spatial_tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), interp_behav_vec, spatial_bin_vector, running_ts);
    [actual_spatial_state_MI(cell_i), ~, spatial_state_occupancy_vector, spatial_state_prob_being_active(cell_i), spatial_state_tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), single_trajectory_vec, cum_spatial_bin_vector, running_ts);
    [actual_distance_MI(cell_i), ~, distance_occupancy_vector, distance_prob_being_active(cell_i), distance_tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), cum_distance, cum_spatial_bin_vector, running_ts);
    [actual_time_MI(cell_i), ~, time_occupancy_vector, time_prob_being_active(cell_i), time_tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), cum_time, cum_temporal_bin_vector, running_ts);
    
    disp([num2str(round(cell_i/numNeurons*100)) '% completed'])
end

save([workingDir filesep 'LMT_cells.mat'], 'spatial_tuning_curve_data', 'spatial_state_tuning_curve_data', 'distance_tuning_curve_data',  'time_tuning_curve_data');

%% Sort data
% Time cells
timeCells_2remove = time_zMI < 2 | spatial_zMI>=2;
sigTime_tuning_curve_data = time_tuning_curve_data;
sigTime_tuning_curve_data(:,timeCells_2remove) = [];
[~, peak_time_index] = max(sigTime_tuning_curve_data,[],1);
[~, sorted_time_index] = sort(peak_time_index);
sorted_time_tuningCurves = sigTime_tuning_curve_data(:,sorted_time_index);

% State-dependent place cells
statePlaceCells_2remove = spatial_state_zMI < 2 | time_zMI>=2 | spatial_zMI>=2;
sigSpaceState_tuning_curve_data = spatial_state_tuning_curve_data;
sigSpaceState_tuning_curve_data(:,statePlaceCells_2remove) = [];
[~, peak_spaceState_index] = max(sigSpaceState_tuning_curve_data,[],1);
[~, sorted_spaceState_index] = sort(peak_spaceState_index);
sorted_spaceState_tuningCurves = sigSpaceState_tuning_curve_data(:,sorted_spaceState_index);

% Place cells
placeCells_2remove = spatial_zMI < 2 | time_zMI>=2;
sigSpace_tuning_curve_data = spatial_tuning_curve_data;
sigSpace_tuning_curve_data(:,placeCells_2remove) = [];
[~, peak_space_index] = max(sigSpace_tuning_curve_data,[],1);
[~, sorted_space_index] = sort(peak_space_index);
sorted_space_tuningCurves = sigSpace_tuning_curve_data(:,sorted_space_index);





