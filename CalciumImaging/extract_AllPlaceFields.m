function extractAllPlaceFields
%EXTRACTALLPLACEFIELDS Summary of this function goes here
%   Detailed explanation goes here

%% Parameters
binSize = 3;
sampling_frequency = 30; % This data set has been sampled at 30 images per second
z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise
min_speed_threshold = 5; % 2 cm.s-1

%% Load ms and behav files in current directory
warning all off
workingFolder = pwd;
load([workingFolder filesep 'ms.mat'])
load([workingFolder filesep 'behav.mat'])

behav_time = behav.time/1000;
ca_time = ms.time/1000;
ca_data = ms.RawTraces;
behav_vec = behav.position;

%% Only keep unique time points
[behav_time, IAbehav, ICbehav]=unique(behav_time);
[ca_time, IAms, ICms]=unique(ca_time);
ca_data = ca_data(IAms,:);
behav_vec = behav_vec(IAbehav,:);

numNeurons = size(ca_data,2);

%% Interpolate behavior
interp_behav_vec(:,1) = interpolate_behavior(behav_vec(:,1), behav_time, ca_time); % in the X dimension
interp_behav_vec(:,2) = interpolate_behavior(behav_vec(:,2), behav_time, ca_time); % in the Y dimension
interp_behav_vec(end,:) = interp_behav_vec(end-1,:);

%% Extract velocity
[velocity] = extract_velocity(interp_behav_vec, ca_time);
running_ts = velocity > min_speed_threshold;

%% Binarize
binarizedData = 0*ca_data;
for cell_i = 1:numNeurons
    binarizedData(:,cell_i) = extract_binary(ca_data(:,cell_i), sampling_frequency, z_threshold);
end

%% Create bins
X_bin_vector = 0:binSize:45+binSize;
X_bin_centers_vector = X_bin_vector + binSize/2;
X_bin_centers_vector(end) = [];

Y_bin_vector = 0:binSize:45+binSize;
Y_bin_centers_vector = Y_bin_vector + binSize/2;
Y_bin_centers_vector(end) = [];

%% Extract place field
for cell_i = 1:numNeurons
    [~, ~, ~, ~, PF] = extract_2D_information(binarizedData(:,cell_i), interp_behav_vec, X_bin_vector, Y_bin_vector, running_ts);
    placeFields(:,:,cell_i) = smoothPlaceField(PF);

end

save([workingFolder filesep 'placeFields.mat'], 'placeFields')
end

