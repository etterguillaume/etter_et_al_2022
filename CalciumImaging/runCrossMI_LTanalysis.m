%% Load data here (ms, behav)
%load('/Users/guillaume/Dropbox (Williams Lab)/Guillaume/analyzed_data/calcium_imaging/J20_246/Linear track/LT7/behav.mat')
%load('/Users/guillaume/Dropbox (Williams Lab)/Guillaume/analyzed_data/calcium_imaging/J20_246/Linear track/LT7/ms.mat')

ca_data = ms.RawTraces;
behav_vec = behav.position(:,1);
%% Make sure time vectors contain unique time values
[behav_time, IAbehav, ICbehav]=unique(behav.time/1000);
[ca_time, IAms, ICms]=unique(ms.time/1000);

ca_data = ca_data(IAms,:);
behav_vec = behav_vec(IAbehav);

%% Extract sampling frequency and other parameters
Fs = round(1/mode(diff(ca_time)));
numFrames = length(ca_time);
numNeurons = size(ms.RawTraces,2);

%% Binarize traces
for cell_i = 1:numNeurons
binarizedData(:,cell_i) = extract_binary(ca_data(:,cell_i), Fs, 2);
end

%% Interpolate behavior
[interp_behav_vec] = interpolate_behavior(behav_vec, behav_time, ca_time);
interp_behav_vec(end) = interp_behav_vec(end-1);

%% Compute velocities
% In this example, we will ignore moments when the mouse is immobile
[ velocity ] = extract_velocity(interp_behav_vec, ca_time);

%% Find periods of immobility
% This will be usefull later to exclude immobility periods from analysis
min_speed_threshold = 5; % 2 cm.s-1
running_ts = velocity > min_speed_threshold;

%% Create an inclusion vector to isolate only a specific behavior
runningData = binarizedData(running_ts,:);

[zMI, MI, p_value] = extractCrossMI(runningData, true);

embedding = run_umap(zMI,'n_components',2);


%% Create graph
adjacency_matrix = MI;
adjacency_matrix(p_value < 0.05) = 0;
adjacency_matrix(isnan(adjacency_matrix)) = 0;
max_MI = max(adjacency_matrix(:));
maxLineWidth = 10;
line_factor = 1/max_MI*maxLineWidth;

figure
scatter(embedding(:,1), embedding(:,2),[],'k','.')
for i = 1:numNeurons
    for j = 1:numNeurons
        if adjacency_matrix(i,j) ~= 0
            line(embedding(i,:),embedding(j,:),'Color','k','LineWidth', adjacency_matrix(i,j)*line_factor)
            drawnow
            pause
        end
    end
end
     
