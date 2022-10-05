%% Load data here
load('/Users/guillaume/Dropbox (Williams Lab)/Guillaume/analyzed_data/calcium_imaging/PV986/LT/MonoTrack/PV986_MT4_20190228/behav.mat')
load('/Users/guillaume/Dropbox (Williams Lab)/Guillaume/analyzed_data/calcium_imaging/PV986/LT/MonoTrack/PV986_MT4_20190228/ms.mat')

%% Keep only unique timestamps and extract
[behav_time, IAbehav, ICbehav]=unique(behav.time/1000);
[ca_time, IAms, ICms]=unique(ms.time/1000);

Fs = round(1/mode(diff(ca_time)));
numFrames = length(ca_time);
numNeurons = size(ms.RawTraces,2);

%% Extract binary activity
for cell_i = 1:numNeurons
binarizedData(:,cell_i) = extract_binary(ms.RawTraces(IAms,cell_i), Fs, 2);
end

%% Interpolate behavior data
[interp_behav_vec] = interpolate_behavior(behav.position(IAbehav,1), behav_time, ca_time);

%% Extract distance travelled
% Not currently used
%dist_travelled = abs(diff(interp_behav_vec));
%dist_travelled(end+1) = 0;
%dist_travelled = cumsum(dist_travelled);

%% Extract velocity
[ velocity ] = extract_velocity(interp_behav_vec, ca_time);

runningThreshold = 2;
isRunning = velocity > runningThreshold;

%% Extract end portions
isRight = interp_behav_vec < 20;
isLeft = interp_behav_vec > 134-20;

%% Extract sound cue
[interp_sound_vec] = interpolate_behavior(behav.optosignal(IAbehav)', behav_time, ca_time);

soundThreshold = 0.9;
isSoundOn = interp_sound_vec > soundThreshold;

%% Extract sound cue onset/offset
diff_sound = diff(isSoundOn);
diff_sound(2:end+1) = diff_sound;
diff_sound(1)=0;
isSoundOnset = diff_sound > 0;
isSoundOffset = diff_sound < 0;    

%% Generate states matrix
states_matrix(:,1) = isLeft; % Is the mouse in the leftt portion of the maze
states_matrix(:,2) = isRight; % Is the mouse in the right portion of the maze
states_matrix(:,3) = isSoundOn; % Sound cue
states_matrix(:,4) = isSoundOnset; % Sound onset
states_matrix(:,5) = isSoundOffset; % Sound offset

states_matrix(~isRunning,:) = 0;

%% MI analysis
[MI, p_values, zMI] = MI_state_analysis(binarizedData, states_matrix);
pThreshold = 0.05;
MI(p_values < pThreshold) = NaN;

%% Fit a GLM to predictors/cell activity
%[b, dev, stats] = glmfit(states_matrix,binarizedData(:,1),'binomial');






