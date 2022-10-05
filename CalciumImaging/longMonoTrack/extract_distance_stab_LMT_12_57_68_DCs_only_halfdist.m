%% Parameters
spatialBinSize = 3;
max_cum_spatial_length = 268;
sampling_frequency = 30; % This data set has been sampled at 30 images per second
z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise
min_speed_threshold = 5; %cm.s-1
lap_detect = 10;

%% Create bins
cum_spatial_bin_vector = 0:spatialBinSize:max_cum_spatial_length+spatialBinSize;
cum_spatial_bin_centers_vector = cum_spatial_bin_vector + spatialBinSize/2;
cum_spatial_bin_centers_vector(end) = [];

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
spatial_pvalues = {};
optosignal = {};

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
        optosignal{sesh} = behav.optosignal;
        load([workingFolder filesep session_folders(currentFolder).name '/LMT_results.mat'])
        spatial_pvalues{sesh} = spatial_MI_pvalue;
        temporal_pvalues{sesh} = time_MI_pvalue;
        distance_pvalues{sesh} = distance_MI_pvalue;
        sesh = sesh +1;
    end
end

%% Preprocess the data
interp_behav_vec = {};
LMT_state = {};
cum_distance = {};
binarizedData = {};
running_ts = {};
Fs = sampling_frequency;

for session_i = 1:length(ca_data)
    % Only keep unique time points
    LMT_state{session_i} = extractLMTState(optosignal{session_i});
    [behav_time{session_i}, IAbehav, ~] = unique(behav_time{session_i});
    [ca_time{session_i}, IAms, ~] = unique(ca_time{session_i});
    ca_data{session_i} = ca_data{session_i}(IAms,:);
    behav_vec{session_i} = behav_vec{session_i}(IAbehav,:);
    LMT_state{session_i} = LMT_state{session_i}(IAbehav);
    numNeurons{session_i} = size(ca_data{session_i},2);
    
    %% Interpolate behavior
    interp_behav_vec{session_i} = interpolate_behavior(behav_vec{session_i}(:,1), behav_time{session_i}, ca_time{session_i}); % in the X dimension
    interp_behav_vec{session_i}(end) = interp_behav_vec{session_i}(end-1);
    
    %% Extract state
    LMT_state{session_i} = interp1(behav_time{session_i}, LMT_state{session_i},ca_time{session_i},'nearest');

    %% Measure elapsed time
    travelled_dist{session_i} = diff(interp_behav_vec{session_i});
    travelled_dist{session_i}(end+1)=0;
    
    cum_distance{session_i}(1) = 0;
    
    for frame_i = 2:length(LMT_state{session_i})
        if LMT_state{session_i}(frame_i) == 1 && interp_behav_vec{session_i}(frame_i) < lap_detect
            cum_distance{session_i}(frame_i) = 0;
        else    
            cum_distance{session_i}(frame_i) = cum_distance{session_i}(frame_i-1)+abs(travelled_dist{session_i}(frame_i));
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

%% Identify place cells
distanceCellsProper_1 = spatial_pvalues{1} > 0.05 & temporal_pvalues{1} > 0.05 & distance_pvalues{1}<0.05;
distanceCellsProper_3 = spatial_pvalues{3} > 0.05 & temporal_pvalues{3} > 0.05 & distance_pvalues{3}<0.05;
distanceCellsProper_4 = spatial_pvalues{4} > 0.05 & temporal_pvalues{4} > 0.05 & distance_pvalues{4}<0.05;

%% Identify cell presents in sessions of interest
cellIndex_12 = [];
cellIndex_35 = [];
cellIndex_46 = [];
for cell_i = 1:size(cell_to_index_map)
    if cell_to_index_map(cell_i,1)>0 && cell_to_index_map(cell_i,2)>0 % Select only cells active from day 3 (place cells) to 5
        if distanceCellsProper_1(cell_to_index_map(cell_i,1)) == 1
            cellIndex_12(end+1,:) = cell_i;
        end
    end
    if cell_to_index_map(cell_i,3)>0 && cell_to_index_map(cell_i,5)>0 % Select only cells active from day 3 (place cells) to 5
        if distanceCellsProper_3(cell_to_index_map(cell_i,3)) == 1
            cellIndex_35(end+1,:) = cell_i;
        end
    end
    if cell_to_index_map(cell_i,4)>0 && cell_to_index_map(cell_i,6)>0 % Select only cells active from day 4 (place cells) to 6
        if distanceCellsProper_4(cell_to_index_map(cell_i,4)) == 1
            cellIndex_46(end+1,:) = cell_i;
        end
    end
end

%% Assess stability across days
%% Day 1-2
ct=1;
distanceFields1=[];
distanceFields2=[];
for cell_i = 1:length(cellIndex_12)
    [~, ~, ~, ~, DF_1] = extract_1D_information(binarizedData{1}(:,cell_to_index_map(cellIndex_12(cell_i),1)), cum_distance{1}, cum_spatial_bin_vector, running_ts{1});
    DF_1 = smoothPlaceField(DF_1);
    DF_1 = DF_1-min(DF_1);
    DF_1 = DF_1./max(DF_1);
    distanceFields1(:,end+1) = DF_1;
    
    [~, ~, ~, ~, PF_2] = extract_1D_information(binarizedData{2}(:,cell_to_index_map(cellIndex_12(cell_i),2)), cum_distance{2}, cum_spatial_bin_vector, running_ts{2});
    PF_2 = smoothPlaceField(PF_2);
    PF_2 = PF_2-min(PF_2);
    PF_2 = PF_2./max(PF_2);
    distanceFields2(:,end+1) = PF_2;

    stability_12(cell_i) = corr(PF_2,DF_1);
    ct = ct+1;   
end

%% Day 3-5
ct=1;
distanceFields3=[];
distanceFields5=[];
for cell_i = 1:length(cellIndex_35)
    [~, ~, ~, ~, PF_3] = extract_1D_information(binarizedData{3}(:,cell_to_index_map(cellIndex_35(cell_i),3)), cum_distance{3}, cum_spatial_bin_vector, running_ts{3});
    PF_3 = smoothPlaceField(PF_3);
    PF_3 = PF_3-min(PF_3);
    PF_3 = PF_3./max(PF_3);
    distanceFields3(:,end+1) = PF_3;
    
    [~, ~, ~, ~, PF_5] = extract_1D_information(binarizedData{5}(:,cell_to_index_map(cellIndex_35(cell_i),5)), cum_distance{5}, cum_spatial_bin_vector, running_ts{5});
    PF_5 = smoothPlaceField(PF_5);
    PF_5 = PF_5-min(PF_5);
    PF_5 = PF_5./max(PF_5);
    distanceFields5(:,end+1) = PF_5;

    stability_35(cell_i) = corr(PF_5,PF_3);
    ct = ct+1;   
end

%% Day 4-6
ct=1;
distanceFields4=[];
distanceFields6=[];
for cell_i = 1:length(cellIndex_46)
    [~, ~, ~, ~, PF_4] = extract_1D_information(binarizedData{4}(:,cell_to_index_map(cellIndex_46(cell_i),4)), cum_distance{4}, cum_spatial_bin_vector, running_ts{4});
    PF_4 = smoothPlaceField(PF_4);
    PF_4 = PF_4-min(PF_4);
    PF_4 = PF_4./max(PF_4);
    distanceFields4(:,end+1) = PF_4;
    
    [~, ~, ~, ~, PF_6] = extract_1D_information(binarizedData{6}(:,cell_to_index_map(cellIndex_46(cell_i),6)), cum_distance{6}, cum_spatial_bin_vector, running_ts{6});
    PF_6 = smoothPlaceField(PF_6);
    PF_6 = PF_6-min(PF_6);
    PF_6 = PF_6./max(PF_6);
    distanceFields6(:,end+1) = PF_6;

    stability_46(cell_i) = corr(PF_6,PF_4);
    ct = ct+1;   
end

save([workingFolder filesep 'LMT_distance_stability_12_35_46_DCs_only_halfdist.mat'],'stability_12', 'stability_35','stability_46', 'distanceFields1','distanceFields2','distanceFields3', 'distanceFields4', 'distanceFields5', 'distanceFields6', 'mouse_number')

%% Plot here
figure
%% 1-2
[~, maxIdx] = max(distanceFields1);
[~, sortedIdx] = sort(maxIdx,'descend');

subplot(3,2,1) 
imagesc(distanceFields1(:,sortedIdx)')
colormap Viridis

subplot(3,2,2) 
imagesc(distanceFields2(:,sortedIdx)')
colormap Viridis

%% 3-5
[~, maxIdx] = max(distanceFields3);
[~, sortedIdx] = sort(maxIdx,'descend');

subplot(3,2,3) 
imagesc(distanceFields3(:,sortedIdx)')
colormap Viridis

subplot(3,2,4) 
imagesc(distanceFields5(:,sortedIdx)')
colormap Viridis
%% 4-6
[~, maxIdx] = max(distanceFields4);
[~, sortedIdx] = sort(maxIdx,'descend');

subplot(3,2,5) 
imagesc(distanceFields4(:,sortedIdx)')
colormap Viridis

subplot(3,2,6) 
imagesc(distanceFields6(:,sortedIdx)')
colormap Viridis