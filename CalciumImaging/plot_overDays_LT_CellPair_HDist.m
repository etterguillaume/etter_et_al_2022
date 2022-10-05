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
    numNeurons=size(ca_data{session_i},2);
    
    %% Interpolate behavior
    interp_behav_vec{session_i} = interpolate_behavior(behav_vec{session_i}(:,1), behav_time{session_i}, ca_time{session_i}); % in the X dimension
    interp_behav_vec{session_i}(end) = interp_behav_vec{session_i}(end-1);

    %% Extract velocity
    velocity{session_i} = extract_velocity(interp_behav_vec{session_i}, ca_time{session_i});
    
    %% Binarize    
    binarizedData{session_i} = 0*ca_data{session_i};
    for cell_i = 1:numNeurons
        binarizedData{session_i}(:,cell_i) = extract_binary(ca_data{session_i}(:,cell_i), sampling_frequency, z_threshold);
    end
    
    %% Keep only running epochs
    running_ts = velocity{session_i} > min_speed_threshold;
    binarizedData{session_i} = binarizedData{session_i}(running_ts,:);
    
end

%% Identify cell presents in sessions of interest
cellIndex = [];
for cell_i = 1:size(cell_to_index_map)
    if sum(cell_to_index_map(cell_i,:)==0) == 0
        cellIndex(end+1,:) = cell_i;
    end
end

numNeurons = length(cellIndex);

%% Pre-allocate memory and create labels
overDaysLabels = {'Day 1','Day 2', 'Day 3', 'Day 4', 'Day 5', 'Day 6', 'Day 7, scrambled stims', 'Day 8, 8Hz stims', 'Day 9', 'Day 10'};

HDist={}
figure
for currentDay = 1:10
    disp(['Current day: ' num2str(currentDay) '/10'])
    data = binarizedData{currentDay}(:,cell_to_index_map(cellIndex,currentDay));
    
    for cell_i = 1:numNeurons
        disp(['Progress: ' num2str(cell_i/numNeurons*100) '%']);
        for cell_j = 1:numNeurons

            if cell_j < cell_i % This is to avoid computing the symmetric part (and save time)
                % First compute marginal probabilities that will not change
                % after shuffling:
                trace_i = data(:,cell_i);
                trace_j = data(:,cell_j);

                HDist{currentDay}(cell_i,cell_j) = extract_HammingDistance(trace_i,trace_j);
                
            end
        end
    end
    
    subplot(1,10,currentDay)
    imagesc(HDist{currentDay})
    ax=gca;
    ax.CLim=[0 1];
    daspect([1 1 1])
    colorbar
    colormap Viridis
    drawnow

end

save([workingFolder filesep 'HDist_to_plot.mat'],'HDist', 'mouse_number')