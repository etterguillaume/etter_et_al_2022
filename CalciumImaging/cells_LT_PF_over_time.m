function results = cells_LT_PF_over_time(sessions2include)
%% This function will extract place fields from linear track sessions for cells that were detected in all included sessions

%% Parameters used
% Bootstrapping
numFrames2Use = 500;
numShuffles = 200;

% Bin size
bin_size = 3;
    % Make sure that your binning vector includes every data point of
    % interp_behav_vec using the min/max function:
    bin_vector = 0:bin_size:134+bin_size; % start : bin_size : end
    bin_centers_vector = bin_vector + bin_size/2;
    bin_centers_vector(end) = [];

%% Load the registration results from CellReg
load cellRegistered.mat
cell2index = cell_registered_struct.cell_to_index_map;

% ADD ERROR CHECKING FOR sessions2include VECTOR
% Order if disordered
sessions2include = sort(sessions2include,'ascend');

if isempty(sessions2include)
    sessions2include = 1:size(cell2index,2); % Include every session
end

%% Define the cells that will be used as being detected in all included sessions
cells2use = sum(cell2index(:,sessions2include) == 0,2) == 0;
% Now we can remove all the unused_cells in cell2index
cell2index = cell2index(cells2use,:);

%% Now let's list all the sessions and sort them in chronological order
folders = dir('*LT*');
for i=1:length(folders)
    folderNames{i} = folders(i).name;
end

folderList = [];
dateList = [];
for i = 1:length(folderNames)
    if contains(folderNames{i},'LT')
        folderList{end+1} = folderNames{i};
        temp = split(folderNames{i},'_');
        dateList(end+1) = str2num(temp{end});
    end
end

[~,sortInd] = sort(dateList,'ascend');
sortedFolderList = folderList(sortInd);

%% Now let's load included sessions one after the other
warning off all % Removes warnings associated with Matlab attempting to load every ms video

results = {};
for session_i = 1:length(sessions2include)
    disp(['Session: ' num2str(session_i) '/' num2str(length(sessions2include))])
    load([pwd filesep sortedFolderList{session_i} filesep 'ms.mat']);
    load([pwd filesep sortedFolderList{session_i} filesep 'behav.mat']);
    
    ca_data = ms.RawTraces(:,cell2index(:,sessions2include(session_i))); % Perfom computations only on included cells
    ca_time = ms.time/1000;
    
    behav_vec = behav.position(:,1);
    behav_time = behav.time/1000;
    
    %% Keep only unique data points
    [behav_time, IAbehav, ICbehav]=unique(behav_time);
    behav_vec = behav_vec(IAbehav);
    [ca_time, IAms, ICms]=unique(ca_time);
    ca_data = ca_data(IAms,:);
    
    %% Binarize data
    sampling_frequency = 30; % This data set has been sampled at 30 images per second
    z_threshold = 2;
    numNeurons = size(ca_data,2);
    
    binarized_data = zeros(size(ca_data));
    for cell_i = 1:numNeurons
        binarized_data(:,cell_i) = extract_binary(ca_data(:,cell_i), sampling_frequency, z_threshold);
    end
    
    %% Interpolate behavior
    [interp_behav_vec] = interpolate_behavior(behav_vec, behav_time, ca_time);
    interp_behav_vec(end) = interp_behav_vec(end-1);
    
    %% Compute locomotor velocity
    [ velocity ] = extract_velocity(interp_behav_vec, ca_time);
    
    %% Identify running directions
    [right_indices] = isolate_direction(interp_behav_vec,'right');
    [left_indices] = isolate_direction(interp_behav_vec,'left');
    
    %% Start building an inclusion vector that only considers running periods
    min_speed_threshold = 5; % 2 cm.s-1
    running_ts = velocity > min_speed_threshold;
    
    %% Create inclusion vectors for right, left, and both trajectories
    both_inclusion_vector = running_ts;
    right_inclusion_vector = right_indices & running_ts;
    left_inclusion_vector = left_indices & running_ts;
    
    %% Error checking if sampling is too low:
    if sum(right_inclusion_vector) < 2*numFrames2Use
        warning('The number of right trajectory running epochs is low. Bootstrap samples will be lowered consequently')
        numFrames2Use = round(sum(right_inclusion_vector)/2);
    elseif sum(left_inclusion_vector) < 2*numFrames2Use
        warning('The number of left trajectory running epochs is low. Bootstrap samples will be lowered consequently')
        numFrames2Use = round(sum(left_inclusion_vector)/2);
    end
    
    %% Now we will create a m x n matrix where m represents the location on the linear track and n the bootstrap sample
    actual_bootstrap_tuning_curve = zeros(length(bin_centers_vector), numShuffles);
    
    for k = 1:numShuffles
        disp(['Bootstrap sample: ' num2str(k) '/' num2str(numShuffles)])
        both_bootstrap_ts = zeros(size(both_inclusion_vector));
        right_bootstrap_ts = zeros(size(both_inclusion_vector));
        left_bootstrap_ts = zeros(size(both_inclusion_vector));
        
        both_included_loc = find(both_inclusion_vector == 1);
        right_included_loc = find(right_inclusion_vector == 1);
        left_included_loc = find(left_inclusion_vector == 1);
        
        for frame_i = 1:numFrames2Use
            both_random_frame = ceil(rand*length(both_included_loc));
            right_random_frame = ceil(rand*length(right_included_loc));
            left_random_frame = ceil(rand*length(left_included_loc));
            
            both_bootstrap_ts(both_included_loc(both_random_frame)) = 1;
            right_bootstrap_ts(right_included_loc(right_random_frame)) = 1;
            left_bootstrap_ts(left_included_loc(left_random_frame)) = 1;
        end
        
        both_bootstrap_ts = logical(both_bootstrap_ts);
        right_bootstrap_ts = logical(right_bootstrap_ts);
        left_bootstrap_ts = logical(left_bootstrap_ts);
        
        %% Extract spatial information and append to the results cell
        for cell_i = 1:numNeurons
            [results{session_i,cell_i}.both_MI(k), results{session_i,cell_i}.both_posterior(:,k), results{session_i,cell_i}.both_occupancy_vector, results{session_i,cell_i}.both_prob_being_active(k), results{session_i,cell_i}.both_tuning_curve(:,k) ] = extract_1D_information(binarized_data(:,cell_i), interp_behav_vec, bin_vector, both_bootstrap_ts);
            [results{session_i,cell_i}.right_MI(k), results{session_i,cell_i}.right_posterior(:,k), results{session_i,cell_i}.right_occupancy_vector, results{session_i,cell_i}.right_prob_being_active(k), results{session_i,cell_i}.right_tuning_curve(:,k) ] = extract_1D_information(binarized_data(:,cell_i), interp_behav_vec, bin_vector, right_bootstrap_ts);
            [results{session_i,cell_i}.left_MI(k), results{session_i,cell_i}.left_posterior(:,k), results{session_i,cell_i}.left_occupancy_vector, results{session_i,cell_i}.left_prob_being_active(k), results{session_i,cell_i}.left_tuning_curve(:,k) ] = extract_1D_information(binarized_data(:,cell_i), interp_behav_vec, bin_vector, left_bootstrap_ts);
        end
    end

end
end
