workingFolder = pwd;

%% Parameters
numShuffles = 1000;
binSize = 3;
sampling_frequency = 30; % This data set has been sampled at 30 images per second
z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise
min_speed_threshold = 5; % 2 cm.s-1

%% Create bins
bin_vector = 0:binSize:134*2+binSize;
bin_centers_vector = bin_vector + binSize/2;
bin_centers_vector(end) = [];

%% Load the data
load([workingFolder filesep 'ms.mat'])
ca_data = ms.RawTraces;
ca_time = ms.time/1000;
load([workingFolder filesep 'behav.mat'])
behav_vec = behav.position;
behav_time = behav.time/1000;

%% Only keep unique time points
[behav_time, IAbehav, ~] = unique(behav_time);
[ca_time, IAms, ~] = unique(ca_time);
ca_data = ca_data(IAms,:);
behav_vec = behav_vec(IAbehav,:);
numNeurons = size(ca_data,2);

%% Interpolate behavior
interp_behav_vec = interpolate_behavior(behav_vec(:,1), behav_time, ca_time); % in the X dimension
interp_behav_vec(end) = interp_behav_vec(end-1);

%% Extract velocity
velocity = extract_velocity(interp_behav_vec, ca_time);

%% Binarize
binarizedData = 0*ca_data;
for cell_i = 1:numNeurons
    binarizedData(:,cell_i) = extract_binary(ca_data(:,cell_i), sampling_frequency, z_threshold);
end

running_ts = velocity > min_speed_threshold;

%% Extract actual MI for each cell
for cell_i = 1:numNeurons
    [actualMI(cell_i), ~, ~, ~, ~] = extract_1D_information(binarizedData(:,cell_i), interp_behav_vec, bin_vector, running_ts);
end

%% Shuffle data
shuffledMI = zeros(numNeurons,numShuffles);

for k = 1:numShuffles
    disp(k)
    random_ts = ceil(rand*length(ca_time));
    shuffled_binarized = zeros(length(ca_time),numNeurons);

    % Permute the trace
    shuffled_binarized(1:random_ts,:) = binarizedData(end-random_ts+1:end,:);
    shuffled_binarized(end-random_ts+1:end,:) = binarizedData(1:random_ts,:);
    
    % Compute tuning curve
    for cell_i = 1:numNeurons
        [shuffledMI(cell_i,k), ~, ~, ~, ~] = extract_1D_information(shuffled_binarized(:,cell_i), interp_behav_vec, bin_vector, running_ts);
    end
end 

%% Compute significance
significantPFs = zeros(numNeurons,1);
for cell_i = 1:numNeurons
    pvalue = sum(shuffledMI(cell_i,:) > actualMI(cell_i),2)/numShuffles; %  p-value, supra-threshold tests
    if pvalue < 0.05
        significantPFs(cell_i) = 1;
    end
end

save([workingFolder filesep 'significantPFs'], 'significantPFs')