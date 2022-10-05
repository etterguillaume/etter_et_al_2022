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

sesh = 1;
for folder_i = 1:length(sorted_dates)
    currentFolder = sorted_idx(folder_i);
    if folder2process(currentFolder)
        load([workingFolder filesep session_folders(currentFolder).name '/ms.mat'])
        ca_data{sesh} = ms.RawTraces;
        ca_time{sesh} = ms.time/1000;
        sesh = sesh +1;
    end
end

%% Preprocess the data
binarizedData = {};
running_ts = {};
Fs = 30;
z_threshold = 2;

for session_i = 1:length(ca_data)
    % Only keep unique time points
    [ca_time{session_i}, IAms, ~] = unique(ca_time{session_i});
    ca_data{session_i} = ca_data{session_i}(IAms,:);
    numNeurons{session_i} = size(ca_data{session_i},2);

    %% Binarize
    binarizedData{session_i} = 0*ca_data{session_i};
    for cell_i = 1:numNeurons{session_i}
        binarizedData{session_i}(:,cell_i) = extract_binary(ca_data{session_i}(:,cell_i), Fs, z_threshold);
        %binarizedData{session_i}(:,cell_i) = extract_binary_MAD(ca_data{session_i}(:,cell_i), Fs);
        %[~,deconvolved_ca,~] = deconvolveCa(ca_data{session_i}(:,cell_i));
        %binarizedData{session_i}(:,cell_i) = deconvolved_ca>0;
    end
    
end

%% Identify cell presents in sessions of interest
cellIndex = [];
for cell_i = 1:size(cell_to_index_map)
    if sum(cell_to_index_map(cell_i,6:9)==0) == 0
        cellIndex(end+1,:) = cell_i;
    end
end

%% Plot cells
figure
for cell_i = 1:length(cellIndex)
    for currentDay = 6:9
        subplot(4,1,currentDay-5)
        plot(ca_time{currentDay},ca_data{currentDay}(:,cell_to_index_map(cellIndex(cell_i),currentDay)),'color','k')
        hold on
        plot(ca_time{currentDay},-binarizedData{currentDay}(:,cell_to_index_map(cellIndex(cell_i),currentDay))+max(ca_data{currentDay}(:,cell_to_index_map(cellIndex(cell_i),currentDay)))+1,'color','red');
        ax(currentDay)=gca;
        ax(currentDay).YLim=[0 6];
        linkaxes(ax,'x','y')
        hold off
    end
    drawnow
    pause
end
