%% Parameters
sampling_frequency = 30; % This data set has been sampled at 30 images per second
z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise
min_speed_threshold = 5; % 2 cm.s-1

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
binarizedData = {};
interp_behav_vec={};
velocity={};
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
    
    %% Keep only running epochs
    running_ts = velocity{session_i} > min_speed_threshold;
    binarizedData{session_i} = binarizedData{session_i}(running_ts,:);
    
end

%% Pre-allocate memory and create labels
eachDay_LT_HDist = {};
eachDayLabels = {'Day 1','Day 2', 'Day 3', 'Day 4', 'Day 5', 'Day 6', 'Day 7, scrambled stims', 'Day 8, 8Hz stims', 'Day 9', 'Day 10'};

for currentDay = 1:10
    disp(['Current day: ' num2str(currentDay) '/10'])
    for cell_i = 1:numNeurons{currentDay}
        for cell_j = 1:numNeurons{currentDay}
            if cell_j > cell_i % to avoid redundancy
                eachDay_LT_HDist{currentDay}(cell_i,cell_j) = extract_unbiasedHammingDistance(binarizedData{currentDay}(:,cell_i), binarizedData{currentDay}(:,cell_j)); 
            else
                eachDay_LT_HDist{currentDay}(cell_i,cell_j) = NaN; 
            end
        end
    end
end

save([workingFolder filesep 'eachDay_LT_CellPair_unbiasedHDist.mat'],'eachDay_LT_HDist', 'eachDayLabels', 'mouse_number')