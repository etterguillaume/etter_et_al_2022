%% Parameters
binSize = 3;
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

%% Assess stability across days
%% Day 3-5
ct=1;
placeFields3 = [];
placeFields5 = [];
for cell_i = 1:length(cellIndex_35)
    [~, ~, ~, ~, PF_3] = extract_1D_information(binarizedData{3}(:,cell_to_index_map(cellIndex_35(cell_i),3)), interp_behav_vec{3}, bin_vector, running_ts{3});
    PF_3 = smoothPlaceField(PF_3);
    PF_3 = PF_3-min(PF_3);
    PF_3 = PF_3./max(PF_3);
    placeFields3(:,end+1) = PF_3;
    
    [~, ~, ~, ~, PF_5] = extract_1D_information(binarizedData{5}(:,cell_to_index_map(cellIndex_35(cell_i),5)), interp_behav_vec{5}, bin_vector, running_ts{5});
    PF_5 = smoothPlaceField(PF_5);
    PF_5 = PF_5-min(PF_5);
    PF_5 = PF_5./max(PF_5);
    placeFields5(:,end+1) = PF_5;

    stability_35(cell_i) = corr(PF_5,PF_3);
    ct = ct+1;   
end

%% Day 4-6
ct=1;
placeFields4 = [];
placeFields6 = [];
for cell_i = 1:length(cellIndex_46)
    [~, ~, ~, ~, PF_4] = extract_1D_information(binarizedData{4}(:,cell_to_index_map(cellIndex_46(cell_i),4)), interp_behav_vec{4}, bin_vector, running_ts{4});
    PF_4 = smoothPlaceField(PF_4);
    PF_4 = PF_4-min(PF_4);
    PF_4 = PF_4./max(PF_4);
    placeFields4(:,end+1) = PF_4;
    
    [~, ~, ~, ~, PF_6] = extract_1D_information(binarizedData{6}(:,cell_to_index_map(cellIndex_46(cell_i),6)), interp_behav_vec{6}, bin_vector, running_ts{6});
    PF_6 = smoothPlaceField(PF_6);
    PF_6 = PF_6-min(PF_6);
    PF_6 = PF_6./max(PF_6);
    placeFields6(:,end+1) = PF_6;

    stability_46(cell_i) = corr(PF_6,PF_4);
    ct = ct+1;   
end

save([workingFolder filesep 'LMT_spatial_stability_35_46.mat'],'stability_35','stability_46', 'mouse_number')

%% Plot here
figure
%% 3-5
[~, maxIdx] = max(placeFields3);
[~, sortedIdx] = sort(maxIdx,'descend');

subplot(2,2,1) 
imagesc(placeFields3(:,sortedIdx)')
colormap Viridis

subplot(2,2,2) 
imagesc(placeFields5(:,sortedIdx)')
colormap Viridis
%% 4-6
[~, maxIdx] = max(placeFields4);
[~, sortedIdx] = sort(maxIdx,'descend');

subplot(2,2,3) 
imagesc(placeFields4(:,sortedIdx)')
colormap Viridis

subplot(2,2,4) 
imagesc(placeFields6(:,sortedIdx)')
colormap Viridis