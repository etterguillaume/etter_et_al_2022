function extract_OF_8Hz_stability

%% Parameters
binSize = 3;
Fs = 30; % This data set has been sampled at 30 images per second
z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise
min_speed_threshold = 2; %cm.s-1

%% Create bins
X_bin_vector = 0:binSize:50+binSize;
X_bin_centers_vector = X_bin_vector + binSize/2;
X_bin_centers_vector(end) = [];

Y_bin_vector = 0:binSize:50+binSize;
Y_bin_centers_vector = Y_bin_vector + binSize/2;
Y_bin_centers_vector(end) = [];

%% Open files
workingDir = pwd;

load([workingDir filesep 'ms.mat'])
load([workingDir filesep 'behav.mat'])

ca_time = ms.time/1000;
ca_data = ms.RawTraces;

behav_vec = behav.position;
behav_time = behav.time/1000;

[behav_time, IAbehav, ~] = unique(behav_time);
[ca_time, IAms, ~] = unique(ca_time);
ca_data = ca_data(IAms,:);
behav_vec = behav_vec(IAbehav,:);
numNeurons = size(ca_data,2);
numFrames = size(ca_data,1);

%% Extract scrambled epochs
[stimulatedEpochs,OpticSigThreshold] = extractOptoStim8Hz(ms,behav);

%% Interpolate behavior
interp_behav_vec = interpolate_behavior(behav_vec, behav_time, ca_time);
interp_behav_vec(end) = interp_behav_vec(end-1);

velocity = extract_velocity(interp_behav_vec, ca_time);
running_ts = velocity > min_speed_threshold;

%% Binarize
binarized_data = 0*ca_data;
for cell_i = 1:numNeurons
    binarized_data(:,cell_i) = extract_binary(ca_data(:,cell_i), Fs, z_threshold);
end

%% Sampling random epochs
sample_BL = zeros(1,numFrames);
sample_STIM = zeros(1,numFrames);
for frame = 1:numFrames
    if running_ts(frame)
        if ~stimulatedEpochs(frame)
            sample_BL(frame) = 1;
        elseif stimulatedEpochs(frame)
            sample_STIM(frame) = 1;
        end
    end
end

sample_BL = logical(sample_BL);
sample_STIM = logical(sample_STIM);

%% Buffer matrices with cell number
placeFields_BL = zeros(length(X_bin_vector)-1, length(Y_bin_vector)-1, numNeurons)*nan;
placeFields_STIM = zeros(length(X_bin_vector)-1, length(Y_bin_vector)-1, numNeurons)*nan;
activity_BL = zeros(1,numNeurons)*nan;
activity_STIM = zeros(1,numNeurons)*nan;
stability = zeros(1,numNeurons)*nan;

%% Create tuning curves
for cell_i = 1:numNeurons
    [~, ~, ~, activity_BL(cell_i), PF_A ] = extract_2D_information(binarized_data(:,cell_i), interp_behav_vec, X_bin_vector, Y_bin_vector, sample_BL);
    PF_A = smoothPlaceField(PF_A);
    PF_A = PF_A-min(PF_A(:));
    placeFields_BL(:,:,cell_i) = PF_A./max(PF_A(:),[],'omitnan');

    [~, ~, ~, activity_STIM(cell_i), PF_B ] = extract_2D_information(binarized_data(:,cell_i), interp_behav_vec, X_bin_vector, Y_bin_vector, sample_STIM);
    PF_B = smoothPlaceField(PF_B);
    PF_B = PF_B-min(PF_B(:));
    placeFields_STIM(:,:,cell_i) = PF_B./max(PF_B(:),[],'omitnan');
    
    stability(cell_i) = corr(PF_A(:),PF_B(:));
end

save([workingDir filesep 'OF_stability_8Hz.mat'],'stability', 'activity_BL','activity_STIM', 'placeFields_BL','placeFields_STIM', 'OpticSigThreshold')

end