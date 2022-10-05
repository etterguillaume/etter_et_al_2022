%% Parameters
numBootstrapSamples = 50;
binSize = 3;
sampling_frequency = 30; % This data set has been sampled at 30 images per second
z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise
min_speed_threshold = 5; %cm.s-1
stability_threshold = 0.3;

%% Create bins
bin_vector = 0:binSize:134+binSize;
bin_centers_vector = bin_vector + binSize/2;
bin_centers_vector(end) = [];

%% Indicate working folder here
workingFolder = pwd;
split_path = strsplit(workingFolder,filesep);
mouse_number = split_path{end-2};

%% Analyze folder structure
session_dates = [];
session_folders = dir(workingFolder);
folder2process = [];
for file_i = 1:length(session_folders)
    split_session = strsplit(session_folders(file_i).name,'_');
    if strcmp(split_session{1},mouse_number)
        session_dates(file_i) = str2num(split_session{end}); % For potential chronological sorting later
        folder2process(file_i) = 1;
    else
        session_dates(file_i) = 0;
        folder2process(file_i) = 0;
    end
    if strcmp(split_session{1},'cellRegistered.mat')
        load([workingFolder filesep session_folders(file_i).name]);
        cell_to_index_map = cell_registered_struct.cell_to_index_map; % Load the identity matrix
    end
end

[sorted_dates, sorted_idx] = sort(session_dates);

ca_data = {};
ca_time = {};
behav_vec = {};
behav_time = {};

sesh = 1;
loadSigFile = 1;
for folder_i = 1:length(sorted_dates)
    currentFolder = sorted_idx(folder_i);
    if folder2process(currentFolder)
        load([workingFolder filesep session_folders(currentFolder).name '/ms.mat'])
        ca_data{sesh} = ms.RawTraces;
        ca_time{sesh} = ms.time/1000;
        load([workingFolder filesep session_folders(currentFolder).name '/behav.mat'])
        behav_vec{sesh} = behav.position;
        behav_time{sesh} = behav.time/1000;
        sesh = sesh +1;
    end
end

%% Preprocess the data
interp_behav_vec = {};
binarizedData = {};
running_ts = {};
Fs = sampling_frequency;

for session_i = 1:length(ca_data)
    % Only keep unique time points
    [behav_time{session_i}, IAbehav, ~] = unique(behav_time{session_i});
    [ca_time{session_i}, IAms, ~] = unique(ca_time{session_i});
    ca_data{session_i} = ca_data{session_i}(IAms,:);
    behav_vec{session_i} = behav_vec{session_i}(IAbehav,:);
    numNeurons{session_i} = size(ca_data{session_i},2);
    
    %% Interpolate behavior
    interp_behav_vec{session_i} = interpolate_behavior(behav_vec{session_i}(:,1), behav_time{session_i}, ca_time{session_i}); % in the X dimension
    interp_behav_vec{session_i}(end) = interp_behav_vec{session_i}(end-1);
    
    %% Extract velocity
    velocity{session_i} = extract_velocity(interp_behav_vec{session_i}, ca_time{session_i});
    
    %% Binarize
    binarizedData{session_i} = 0*ca_data{session_i};
    for cell_i = 1:numNeurons{session_i}
        binarizedData{session_i}(:,cell_i) = extract_binary(ca_data{session_i}(:,cell_i), sampling_frequency, z_threshold);
    end
    running_ts{session_i} = velocity{session_i} > min_speed_threshold;
    
end

%% Identify cell presents in sessions of interest
cellIndex = [];
for cell_i = 1:size(cell_to_index_map)
    if  sum(cell_to_index_map(cell_i,5:10)==0) == 0 % Select only cells active from day 5 to 10
        [~, ~, ~, ~, PF_5] = extract_1D_information(binarizedData{5}(:,cell_to_index_map(cell_i,5)), interp_behav_vec{5}, bin_vector, running_ts{5});
        PF_5 = smoothPlaceField(PF_5);
        PF_5 = PF_5-min(PF_5);
        PF_5 = PF_5./max(PF_5);
        
        [~, ~, ~, ~, PF_10] = extract_1D_information(binarizedData{10}(:,cell_to_index_map(cell_i,10)), interp_behav_vec{10}, bin_vector, running_ts{10});
        PF_10 = smoothPlaceField(PF_10);
        PF_10 = PF_10-min(PF_10);
        PF_10 = PF_10./max(PF_10);
        
        stability = corr(PF_5,PF_10);
        if stability > stability_threshold
            cellIndex(end+1,:) = cell_i;
        end
    end
end

%% Training parameters
training_set_creation_method = 'random'; % 'odd', odd timestamps; 'first_portion', first portion of the recording; 3, 'random' random frames
training_set_portion = 0.5; % Portion of the recording used to train the decoder for method 2 and 3

for shuffle_i = 1:numBootstrapSamples
    training_ts = create_training_set(ca_time{1}, training_set_creation_method, training_set_portion);
    training_ts(running_ts{1} == 0) = 0; % Exclude periods of immobility from the traing set
    
    %% Training decoder
    for cell_i = 1:length(cellIndex)
        [~, ~, occupancy_vector, prob_being_active(cell_i), tuning_curve_data(:,cell_i) ] = extract_1D_information(binarizedData{1}(:,cell_i), interp_behav_vec{1}, bin_vector, training_ts);
    end
    
    %% Decode position
    % First, let us establish the timestamps used for decoding.
    decoding_ts = ~training_ts; % Training timestamps are excluded
    decoding_ts(running_ts{1} == 0) = 0; % Periods of immobility are excluded
    
    % Minimal a priori (use to remove experimental a priori)
    occupancy_vector = occupancy_vector./occupancy_vector*(1/length(occupancy_vector));
    
    % Establish which cells are going to be used in the decoding process
    cell_used = logical(ones(size(tuning_curve_data,2),1)); % Let us use every cell for now
    
    mean_decoding_error = {};
    decoding_agreement = {};
    
    for day = 5:10
        [decoded_probabilities] = bayesian_decode1D(binarizedData{day}, occupancy_vector, prob_being_active, tuning_curve_data, cell_used);
        
        [max_decoded_prob, decoded_bin] = max(decoded_probabilities,[],1);
        decoded_position = bin_centers_vector(decoded_bin);
        
        % Before looking at the error rate, we must first bin the actual data using the same bin vector used by
        % the decoder
        actual_bin = nan*interp_behav_vec{day};
        actual_position = nan*interp_behav_vec{day};
        for bin_i = 1:length(bin_vector)-1
            position_idx = find(interp_behav_vec{day}>bin_vector(bin_i) & interp_behav_vec{day} < bin_vector(bin_i+1));
            actual_bin(position_idx) = bin_i;
            actual_position(position_idx) = bin_centers_vector(bin_i);
        end
        
        decoded_bin(~decoding_ts) = nan;
        decoded_position(~decoding_ts) = nan;
        decoded_probabilities(:,~decoding_ts) = nan;
        
        actual_bin(~decoding_ts) = nan;
        actual_position(~decoding_ts) = nan;
        actual_bin = actual_bin';
        actual_position =  actual_position';
        
        %% Compute decoding agreement
        decoding_agreement_vector = double(decoded_bin == actual_bin);
        decoding_agreement_vector(isnan(decoded_bin)) = nan;
        decoding_agreement_vector(isnan(actual_bin)) = nan;
        decoding_agreement_vector(isnan(decoding_agreement_vector)) = [];
        decoding_agreement{day}(shuffle_i) = sum(decoding_agreement_vector)./length(decoding_agreement_vector);
        
        %% Compute decoding error
        decoding_error = actual_position - decoded_position;
        mean_decoding_error{day}(shuffle_i) = mean(abs(decoding_error), 'omitnan');
    end
end
