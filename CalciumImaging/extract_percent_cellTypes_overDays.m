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

%% Process data
sesh = 1;
loadSigFile = 1;
for folder_i = 1:length(sorted_dates)
    currentFolder = sorted_idx(folder_i);
    if folder2process(currentFolder)
        load([workingFolder filesep session_folders(currentFolder).name '/LMT_results.mat'])
        
        placeCellsProper = spatial_MI_pvalue < 0.05 & time_pvalue > 0.05 & pooled_distance_pvalue>0.05;
        timeCellsProper = time_MI_pvalue < 0.05 & spatial_pvalue>0.05 & pooled_distance_pvalue>0.05;
        distanceCellsProper = distance_MI_pvalue < 0.05 & spatial_pvalue > 0.05 & pooled_time_pvalue > 0.05;

        
        percentTimeCells = 
        
        
        sesh = sesh +1;
    end
end