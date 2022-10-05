binSize = 5;
sampling_frequency = 30; % This data set has been sampled at 30 images per second
z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise
min_speed_threshold = 5; %cm.s-1

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
        if loadSigFile % Load file containing list of significant place cells
            load([workingFolder filesep session_folders(currentFolder).name '/significantPFs.mat']);
            loadSigFile = 0;
        end
        sesh = sesh +1;
    end
end

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
    
    %% Extracting right trajectories
    %[direction_indices] = isolate_direction(interp_behav_vec{session_i},'right');
    
    %% Extract velocity
    velocity = extract_velocity(interp_behav_vec{session_i}, ca_time{session_i});
    
    %% Binarize
    binarizedData{session_i} = 0*ca_data{session_i};
    for cell_i = 1:numNeurons{session_i}
        binarizedData{session_i}(:,cell_i) = extract_binary(ca_data{session_i}(:,cell_i), sampling_frequency, z_threshold);
    end
    
    %inclusion_vector{session_i} = velocity > min_speed_threshold & direction_indices == 1;
    inclusion_vector{session_i} = velocity > min_speed_threshold;
    
end

%% Identify cell presents in sessions of interest
cellIndex = [];
for cell_i = 1:size(cell_to_index_map)
    if cell_to_index_map(cell_i,1) > 0
        cellIndex(end+1,:) = cell_i;
    end
end

%% Plot place field each day
numShuffles = 1000;
overDays_LT_MinMaxFields={};

figure
for currentDay = 1:10
    for cell_i = 1:length(cellIndex)
        if cell_to_index_map(cellIndex(cell_i),currentDay) > 0
            [~, ~, ~, ~, PF] = extract_1D_information(binarizedData{currentDay}(:,cell_to_index_map(cellIndex(cell_i),currentDay)), interp_behav_vec{currentDay}, bin_vector, inclusion_vector{currentDay});
            actualPF = smoothPlaceField(PF);
            
            %% PC significance
            
            shuffled_MI = zeros(length(bin_centers_vector),numShuffles); % Initialize matrix value
            for shuffle_i = 1:numShuffles
                random_ts = ceil(rand*length(binarizedData{currentDay}));
                shuffledBinarized = zeros(length(binarizedData{currentDay}),1);
                shuffledBinarized(1:random_ts) = binarizedData{currentDay}(end-random_ts+1:end,cell_to_index_map(cellIndex(cell_i),currentDay));
                shuffledBinarized(random_ts+1:end) = binarizedData{currentDay}(1:end-random_ts,cell_to_index_map(cellIndex(cell_i),currentDay));
                
                [~, ~, ~, ~, PF] = extract_1D_information(shuffledBinarized, interp_behav_vec{currentDay}, bin_vector, inclusion_vector{currentDay});
                shuffledPF(:,shuffle_i) = smoothPlaceField(PF);
            end
            
            pvalue = sum(shuffledPF > actualPF,2)/numShuffles; %  p-value, supra-threshold tests
            
            PF = actualPF;
            PF = PF-min(PF);
            PF = PF./max(PF);
            PF(pvalue > 0.05) = NaN;
            overDays_LT_MinMaxFields{currentDay}(:,cell_i) = PF;
        else
            overDays_LT_MinMaxFields{currentDay}(:,cell_i) = bin_centers_vector*NaN;
        end
    end
    
    %% Plot
    data = overDays_LT_MinMaxFields{currentDay};
    if currentDay == 1
        [~, maxIdx] = max(data);
        [~, sortedIdx] = sort(maxIdx,'descend');
    end
    subplot(1,10,currentDay)
    imagesc(data(:,sortedIdx)');
    %daspect([1 1 1])
    colormap 'Viridis'
end