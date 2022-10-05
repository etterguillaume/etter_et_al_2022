%% Parameters
temporalBinSize = 1;
max_cum_temporal_length = 50;
sampling_frequency = 30; % This data set has been sampled at 30 images per second
z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise
min_speed_threshold = 5; %cm.s-1
lap_detect = 10;

%% Create bins
cum_temporal_bin_vector = 0:temporalBinSize:max_cum_temporal_length+temporalBinSize;
cum_temporal_bin_centers_vector = cum_temporal_bin_vector + temporalBinSize/2;
cum_temporal_bin_centers_vector(end) = [];

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
cum_time = {};
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
    
    %% Measure elapsed time
    elapsed_time = diff(ca_time{session_i});
    elapsed_time(end+1) = 0;    
    
    cum_time{session_i}(1) = 0;
    
    for frame_i = 2:length(interp_behav_vec{session_i})
        if interp_behav_vec{session_i}(frame_i) < lap_detect | interp_behav_vec{session_i}(frame_i) > 134-lap_detect
            cum_time{session_i}(frame_i) = 0;
        else    
        cum_time{session_i}(frame_i) = cum_time{session_i}(frame_i-1)+elapsed_time(frame_i);
        end 
    end
    
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
cellIndex_35 = [];
cellIndex_46 = [];
for cell_i = 1:size(cell_to_index_map)
    if cell_to_index_map(cell_i,3)>0 && cell_to_index_map(cell_i,5)>0 % Select only cells active from day 5 to 10
        cellIndex_35(end+1,:) = cell_i;
    end
    if cell_to_index_map(cell_i,4)>0 && cell_to_index_map(cell_i,6)>0 % Select only cells active from day 5 to 10
        cellIndex_46(end+1,:) = cell_i;
    end
end

%% Pre-allocate memory and create labels
%% Day 3-5
ct=1;
for cell_i = 1:length(cellIndex_35)
    [~, ~, ~, ~, TF_3] = extract_1D_information(binarizedData{3}(:,cell_to_index_map(cellIndex_35(cell_i),3)), cum_time{3}, cum_temporal_bin_vector, running_ts{3});
    TF_3 = smoothPlaceField(TF_3);
    TF_3 = TF_3-min(TF_3);
    TF_3 = TF_3./max(TF_3);
    
    [~, ~, ~, ~, TF_5] = extract_1D_information(binarizedData{5}(:,cell_to_index_map(cellIndex_35(cell_i),5)), cum_time{5}, cum_temporal_bin_vector, running_ts{5});
    TF_5 = smoothPlaceField(TF_5);
    TF_5 = TF_5-min(TF_5);
    TF_5 = TF_5./max(TF_5);

    stability_35(ct) = corr(TF_3,TF_5);
    ct = ct+1;
end

%% Day 4-6
ct=1;
for cell_i = 1:length(cellIndex_46)
    [~, ~, ~, ~, TF_4] = extract_1D_information(binarizedData{4}(:,cell_to_index_map(cellIndex_46(cell_i),4)), cum_time{4}, cum_temporal_bin_vector, running_ts{4});
    TF_4 = smoothPlaceField(TF_4);
    TF_4 = TF_4-min(TF_4);
    TF_4 = TF_4./max(TF_4);
    
    [~, ~, ~, ~, TF_6] = extract_1D_information(binarizedData{6}(:,cell_to_index_map(cellIndex_46(cell_i),6)), cum_time{6}, cum_temporal_bin_vector, running_ts{6});
    TF_6 = smoothPlaceField(TF_6);
    TF_6 = TF_6-min(TF_6);
    TF_6 = TF_6./max(TF_6);

    stability_46(ct) = corr(TF_4,TF_6);
    ct = ct+1;
end

save([workingFolder filesep 'LMT_temporal_stability_35_46.mat'],'stability_35','stability_46', 'mouse_number')

%% Plot here
%figure
%[~, maxIdx] = max(stableTF{5});
%[~, sortedIdx] = sort(maxIdx,'descend');

%for currentDay = 5:10
%   subplot(1,6,currentDay-4)
%   imagesc(stableTF{currentDay}(:,sortedIdx)')
%   colormap Viridis
%end