ca_data = ms.RawTraces;
behav_vec = behav.speed;
behav_time = behav.time/1000;
ca_time = ms.time/1000;

%% Binarize calcium trace
sampling_frequency = 30; % This data set has been sampled at 30 images per second
z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise

binarized_data = zeros(size(ca_data));
for cell_i = 1:size(ca_data,2)
    binarized_data(:,cell_i) = extract_binary(ca_data(:,cell_i), sampling_frequency, z_threshold);
end

inclusion_vector = logical(ones(length(ca_time),1));

%% Interpolate behavior
% In most cases, behavior data has to be interpolated to match neural temporal
% activity assuming it has a lower sampling frequency than neural activity
[interp_behav_vec] = interpolate_behavior(behav_vec, behav_time, ca_time);

%% Compute occupancy and joint probabilities
bin_vector = [0 0.5 1 1.5 2 3 4 6 8 12 16 24 32]; % start : bin_size : end
bin_size = bin_vector(2) - bin_vector(1);
bin_centers_vector = bin_vector + bin_size/2;
bin_centers_vector(end) = [];

for cell_i = 1:size(binarized_data,2)
    [KL_divergence(cell_i), PDF(:,cell_i), occupancy_vector, prob_being_active(cell_i), tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), velocity, ca_time, bin_vector, inclusion_vector);
end

%% Plot the tunning curves
[~,max_index] = max(tuning_curve_data,[],1);
[~,sorted_index] = sort(max_index);
sorted_tuning_curve_data = tuning_curve_data(:,sorted_index);

figure
imagesc(bin_centers_vector,1:size(ca_data,2),sorted_tuning_curve_data')
daspect([1 1 1])
title 'Neuronal tuning curves'
xlabel 'Speed (cm.s-1)'
ylabel 'Cell ID'