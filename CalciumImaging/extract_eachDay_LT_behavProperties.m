%% Parameters
end_portion = 0.1;
trackLength = 134;

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

behav_vec = {};
behav_time = {};

sesh = 1;
loadSigFile = 1;
for folder_i = 1:length(sorted_dates)
    currentFolder = sorted_idx(folder_i);
    if folder2process(currentFolder)
        load([workingFolder filesep session_folders(currentFolder).name '/behav.mat'])
        behav_vec{sesh} = behav.position;
        behav_time{sesh} = behav.time/1000;
        sesh = sesh +1;
    end
end

%% Preprocess the data

for session_i = 1:length(behav_vec)
    % Only keep unique time points
    [behav_time{session_i}, IAbehav, ~] = unique(behav_time{session_i});
    behav_vec{session_i} = behav_vec{session_i}(IAbehav,:);
    
    %% Count laps
    laps{session_i} = 0;
    
    for time_i = 1:length(behav_time{session_i})-1
            if behav_vec{session_i}(time_i) <= end_portion*trackLength & behav_vec{session_i}(time_i+1) > end_portion*trackLength
                laps{session_i} = laps{session_i}+1;
            end
    end
    %% Extract velocity
    velocity{session_i} = extract_velocity(behav_vec{session_i}, behav_time{session_i});
    histSpeed{session_i} = hist(velocity{session_i},50)/length(velocity{session_i});
    totalDistance{session_i} = sum(abs(diff(behav_vec{session_i}(:,1))))/100; % in meters    
end

%% Pre-allocate memory and create labels
eachDayLabels = {'Day 1','Day 2', 'Day 3', 'Day 4', 'Day 5', 'Day 6', 'Day 7, scrambled stims', 'Day 8, 8Hz stims', 'Day 9', 'Day 10'};


save([workingFolder filesep 'eachDay_LT_behavProperties.mat'],'velocity', 'histSpeed', 'totalDistance', 'laps', 'eachDayLabels', 'mouse_number')