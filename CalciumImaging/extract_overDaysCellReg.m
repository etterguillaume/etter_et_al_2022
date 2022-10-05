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
    if strcmp(split_session{1},'cellRegistered.mat')
        load([workingFolder filesep session_folders(file_i).name]);
        cell_to_index_map = cell_registered_struct.cell_to_index_map; % Load the identity matrix
    end
end

%% Identify cell presents in sessions of interest
numCellsOverDays = zeros(size(cell_to_index_map,2),size(cell_to_index_map,2));
numCellsPerSessions = histcounts(sum(cell_to_index_map>0,2),[0.5:1:10.5]);
for day_i = 1:size(cell_to_index_map,2)
    for day_j = 1:size(cell_to_index_map,2)
        numCellsOverDays(day_i,day_j) = sum(cell_to_index_map(:,day_i) > 0 & cell_to_index_map(:,day_j) > 0);
    end
end

save([workingFolder filesep 'overDays_LT_CellReg.mat'],'numCellsPerSessions', 'numCellsOverDays', 'mouse_number')