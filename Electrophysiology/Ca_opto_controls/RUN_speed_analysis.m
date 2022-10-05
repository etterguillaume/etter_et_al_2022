%% Load the data:
load('/speed_GRIN_controls.mat') 

%% First, extract velocities from DLC extracted position data:
PV1100_baseline_velocity = extract_DLC_velocity(PV1100_baseline);
PV1100_withscope_velocity = extract_DLC_velocity(PV1100_withscope);
PV1102_baseline_velocity = extract_DLC_velocity(PV1102_baseline);
PV1102_withscope_velocity = extract_DLC_velocity(PV1102_withscope);
PV1103_baseline_velocity = extract_DLC_velocity(PV1103_baseline);
PV1103_withscope_velocity = extract_DLC_velocity(PV1103_withscope);

%% Now bin the data using MATLAB hist function:
bin_size = 0.5;
speed_vec = 0:bin_size:14; %Arbitrary speed vector used to compare velocities
speed_vec = speed_vec + bin_size/2; % Convert to bin centers
speed_vec(end) = [];

binned_PV1100_baseline_velocity = hist(PV1100_baseline_velocity,speed_vec)/30;
binned_PV1100_withscope_velocity = hist(PV1100_withscope_velocity,speed_vec)/30;
binned_PV1102_baseline_velocity = hist(PV1102_baseline_velocity,speed_vec)/30;
binned_PV1102_withscope_velocity = hist(PV1102_withscope_velocity,speed_vec)/30;
binned_PV1103_baseline_velocity = hist(PV1103_baseline_velocity,speed_vec)/30;
binned_PV1103_withscope_velocity = hist(PV1103_withscope_velocity,speed_vec)/30;

%% Averages
PV1100_baseline_mean = mean(PV1100_baseline_velocity);
PV1100_withscope_mean = mean(PV1100_withscope_velocity);
PV1102_baseline_mean = mean(PV1102_baseline_velocity);
PV1102_withscope_mean = mean(PV1102_withscope_velocity);
PV1103_baseline_mean = mean(PV1103_baseline_velocity);
PV1103_withscope_mean = mean(PV1103_withscope_velocity);