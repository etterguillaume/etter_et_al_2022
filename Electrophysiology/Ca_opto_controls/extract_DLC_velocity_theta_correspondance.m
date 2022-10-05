function [theta_per_velocity] = extract_DLC_velocity_theta_correspondance(data,velocity)
%EXTRACT_DLC_VELOCITY_THETA_CORRESPONDANCE Summary of this function goes here
%   Detailed explanation goes here

binSize = 4;
speed_bins = 0:binSize:16;
speed_bins_centers = speed_bins+binSize/2;

Fs_ephys = 1000; % Assumed sampling rate
dt_ephys = 1/Fs_ephys;
data_time = dt_ephys;
for i = 2:length(data)
    data_time(i) = data_time(i-1)+dt_ephys;
end

Fs_behav = 30;
dt_behav = 1/Fs_behav;
behav_time = 0;
for i = 2:length(velocity)
    behav_time(i) = behav_time(i-1)+dt_behav;
end

%% Interpolate velocity data
interp_velocity = interp1(behav_time,velocity,data_time);

%% Extract theta power
% Filter parameters
N = 1024; % order
f_n = (1/dt_ephys)/2; % Nyquist frequency
theta_band = [4 12];
theta_band = theta_band./f_n;

[theta_B] = fir1(N,theta_band);

theta = filtfilt(theta_B,1,data);
theta_power = abs(hilbert(theta));
theta_power = rescale(theta_power);

theta_per_velocity=[];
for bin_i = 1:length(speed_bins)-1
    speed_idx = find(interp_velocity>= speed_bins(bin_i) & interp_velocity < speed_bins(bin_i+1));
    theta_per_velocity(end+1) = median(theta_power(speed_idx));
end

save('speed_theta_correspondance.mat','theta_per_velocity','speed_bins_centers');

end

