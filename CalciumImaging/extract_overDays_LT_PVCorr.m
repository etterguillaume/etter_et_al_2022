%% Parameters
bootstrapSamples = 30;
binSize = 3;
sampling_frequency = 30; % This data set has been sampled at 30 images per second
z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise
min_speed_threshold = 5; %cm.s-1

%% Create bins
bin_vector = 0:binSize:134*2+binSize;

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
        if loadSigFile % Load file containing list of significant place cells
           load([workingFolder filesep session_folders(currentFolder).name '/significantPFs.mat']);
           loadSigFile = 0;
        end
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
    
    %% Create a single trajectory
    z_delta_x = diff(interp_behav_vec{session_i});
    z_delta_x(end+1) = 0;
    z_delta_x(isnan(z_delta_x)) = 0;
    z_delta_x = zscore(z_delta_x);
    
    right_trajectories = interp_behav_vec{session_i};
    right_trajectories(z_delta_x<0.2) = NaN;
    
    left_trajectories = interp_behav_vec{session_i};
    left_trajectories(z_delta_x>-0.2) = NaN;
    
    left_trajectories = 2*max(interp_behav_vec{session_i}) - left_trajectories;
    
    single_trajectory_vec{session_i} = right_trajectories;
    single_trajectory_vec{session_i}(~isnan(left_trajectories)) = left_trajectories(~isnan(left_trajectories));

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
    if sum(cell_to_index_map(cell_i,:)==0) == 0
        cellIndex(end+1,:) = cell_i;
    end
end

%% Pre-allocate memory and create labels
overDays_LT_PVCorr = zeros(10,10);
overDaysLabels = {'Day 1','Day 2', 'Day 3', 'Day 4', 'Day 5', 'Day 6', 'Day 7, scrambled stims', 'Day 8, 8Hz stims', 'Day 9', 'Day 10'};

for day_i = 1:10
    disp(['Current day: ' num2str(day_i) '/10'])
    PV1 = [];
    for cell_i = 1:length(cellIndex)
        [~, ~, ~, ~, PF] = extract_1D_information(binarizedData{day_i}(:,cell_to_index_map(cellIndex(cell_i),day_i)), interp_behav_vec{day_i}, bin_vector, running_ts{day_i});
        PF = smoothPlaceField(PF);
        PV1 = [PV1; PF];
    end
    
    for day_j = 1:10
        PV2 = [];
        for cell_i = 1:length(cellIndex)
            [~, ~, ~, ~, PF] = extract_1D_information(binarizedData{day_j}(:,cell_to_index_map(cellIndex(cell_i),day_j)), interp_behav_vec{day_j}, bin_vector, running_ts{day_j});
            PF = smoothPlaceField(PF);
            PV2 = [PV2; PF];
        end
        overDays_LT_PVCorr(day_i,day_j) = corr(PV1,PV2);
    end
end

save([workingFolder filesep 'overDays_LT_PVCorr.mat'],'overDays_LT_PVCorr', 'overDaysLabels', 'mouse_number')