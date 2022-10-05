function [ out ] = analyzeStimSpeed(data, behav_res )
%GE_SPATIALGAMMA Summary of this function goes here
%   Detailed explanation goes here

%% TODO
%Add hilbert analysis for inst. power?

%% Selecting channel
channel=input('Please enter the channel to be analyzed: ');
stim_duration = input('Please indicate the duration of stimulations (eg 2s): ');
stim_frequency = input('Please indicate the frequenchy (Hz) of stimulations (enter 0 if no specific frequency): ');

min_speed = 0;

%% Parameters
display('Extracting data.....');
Fs = data.LFP.SamplingRate;
dt = 1/Fs;

stim_duration = round(stim_duration * Fs);

%% Chronux parameters
params.Fs = Fs;            % sampling frequency
params.tapers = [1 1];    % taper parameter TODO double check p 12!!! TW K default[1 1]
params.pad = 4;           %
params.fpass = [1 14];  % frequency range ; default [3 55]
params.err = 0;           % 0 for no error calculated
params.trialave = 0;      % when true average the trials

%% Extracting event times
display('Detecting TTLs.....');
events=data.Events.TimeStamps{1};
events_int=diff(events);
events_freq=1./events_int;

stims(1)= round(events(1) * Fs);

%% Detecting stim start time
for i=1:length(events_freq)-1;
    if events_freq(i) < 1 && events_freq(i+1)>=1.75; % Detection algorithm: detects sudden increases in TTL frequency ie putatively stimulated epochs
        stims(end+1) = round(events(i+1) * Fs);
    end
end

%% Removing stims that overlap with recording segment edges
if stims(end)+ stim_duration >= round(data.Spread(2) * Fs)
    stims(end)=[];
end

if stims(1) - stim_duration <= round(data.Spread(1) * Fs)
    stims(1)=[];
end

% Filter parameters
N = 1000; % order
f_n = (1/dt)/2; % Nyquist frequency
theta_band = [4 12];

theta_band = theta_band./f_n;

[theta_B] = fir1(N,theta_band);

LFP_trace = data.LFP.Data(:,channel);
LFP_time = (1:length(LFP_trace))*dt;

%% Getting stimulated indices
stim_vector = 0*LFP_trace;
for stim_i = 1:length(stims)
    stim_vector(stims(stim_i):stims(stim_i)+stim_duration) = 1;
end

%position = behav_res.pD.final_body_positions; % Use body positions
behav_time = behav_res.frame_times;
speed = behav_res.smoothed_speed;

speed(speed<min_speed) = NaN;

display('Processing position data.....');
nan_vector = isnan(speed);
vec_length = 1:size(speed,1);

speed(nan_vector,1) = interp1(vec_length(~nan_vector), speed(~nan_vector,1), vec_length(nan_vector));

%% Filtering frequency bands
display('Filtering and extracting power values.....');
theta = filtfilt(theta_B,1,LFP_trace);

%% Hilbert transforms
theta_power = abs(hilbert(theta));

%% Interpolate mouse position to LFP sampling rate
display('Interpolating data.....');
interpolated_speed=interp1(behav_time,speed,LFP_time);

%% Analyse every BL/STIM segment
for stim_i = 1:length(stims)
    % Stim triggered averages (STA)
    out.speed_STA(:,stim_i) = interpolated_speed(stims(stim_i)-stim_duration:stims(stim_i)+stim_duration);
    
    % Mean speed
    out.BL_mean_speed(stim_i) = mean(interpolated_speed(stims(stim_i)-stim_duration:stims(stim_i)),'omitnan');
    out.STIM_mean_speed(stim_i) = mean(interpolated_speed(stims(stim_i):stims(stim_i)+stim_duration),'omitnan');
    
    % Coefficient of variation
    out.BL_speed_CV(stim_i) = std(interpolated_speed(stims(stim_i)-stim_duration:stims(stim_i)),'omitnan')./mean(interpolated_speed(stims(stim_i)-stim_duration:stims(stim_i)),'omitnan');
    out.STIM_speed_CV(stim_i) = std(interpolated_speed(stims(stim_i):stims(stim_i)+stim_duration),'omitnan')./mean(interpolated_speed(stims(stim_i):stims(stim_i)+stim_duration),'omitnan');
    
    %% FFT processing
    % Compute spectra
    [out.BL_spectrum(stim_i,:),out.f]= mtspectrumc(LFP_trace(stims(stim_i)-stim_duration:stims(stim_i)),params);
    [out.STIM_spectrum(stim_i,:),out.f]= mtspectrumc(LFP_trace(stims(stim_i):stims(stim_i)+stim_duration),params);
    
    % Compute peak power and frequency values
    [out.BL_peak_power(stim_i),index] = max(out.BL_spectrum(stim_i,:));
    out.BL_peak_frequency(stim_i) = out.f(index);
    
    [out.STIM_peak_power(stim_i),index] = max(out.STIM_spectrum(stim_i,:));
    out.STIM_peak_frequency(stim_i) = out.f(index);
    
    % Compute theta CV
    out.BL_theta_CV(stim_i) = std(theta_power(stims(stim_i)-stim_duration:stims(stim_i)),'omitnan')./mean(theta_power(stims(stim_i)-stim_duration:stims(stim_i)),'omitnan');
    out.STIM_theta_CV(stim_i) = std(theta_power(stims(stim_i):stims(stim_i)+stim_duration),'omitnan')./mean(theta_power(stims(stim_i):stims(stim_i)+stim_duration),'omitnan');
    
    % Entrainment fidelity
    if stim_frequency ~= 0
    out.BL_entrainment_fidelity(stim_i) = sum(out.BL_spectrum(stim_i, out.f > stim_frequency - 0.5 & out.f < stim_frequency + 0.5))./sum(out.BL_spectrum(stim_i, out.f > stim_frequency - 4 & out.f < stim_frequency + 12));
    out.STIM_entrainment_fidelity(stim_i) = sum(out.STIM_spectrum(stim_i, out.f > stim_frequency - 0.5 & out.f < stim_frequency + 0.5))./sum(out.STIM_spectrum(stim_i, out.f > stim_frequency - 4 & out.f < stim_frequency + 12));
    end
    
end

%% Compute STA for speed values
out.mean_speed_STA = mean(out.speed_STA,2);
out.SEM_speed_STA = std(out.speed_STA,[],2)./sqrt(length(stims));

%% Plot bar graphs of BL vs STIM speed values
out.BL_speed_values = interpolated_speed(stim_vector == 0);
out.STIM_speed_values = interpolated_speed(stim_vector == 1);

%% Compute coefficient of variation of BL vs STIM states
out.total_BL_speed_CV = std(out.BL_speed_values,'omitnan')./mean(out.BL_speed_values,'omitnan');
out.total_STIM_speed_CV = std(out.STIM_speed_values,'omitnan')./mean(out.STIM_speed_values,'omitnan');

%% Compute correlation between speed and frequency for BL and STIM epochs
[~, BL_speed_sorting_idx] = sort(out.BL_mean_speed);
[~, STIM_speed_sorting_idx] = sort(out.STIM_mean_speed);

figure()
subplot(3,1,1)
imagesc(1:59,out.f,out.BL_spectrum(BL_speed_sorting_idx,:)')
colormap('Viridis')
subplot(3,1,2)
plot(out.BL_mean_speed(BL_speed_sorting_idx))
subplot(3,1,3)
scatter(out.BL_mean_speed,out.BL_peak_power)

figure
plot(out.f,mean(out.BL_spectrum));
ax=gca;
ax.XLim=[2 14];
ax.YLim=[0 1.5e-9];


figure()
subplot(3,1,1)
imagesc(1:59,out.f,out.STIM_spectrum(STIM_speed_sorting_idx,:)')
colormap('Viridis')
subplot(3,1,2)
plot(out.STIM_mean_speed(STIM_speed_sorting_idx))
subplot(3,1,3)
scatter(out.STIM_mean_speed,out.STIM_peak_power)

figure
plot(out.f,mean(out.STIM_spectrum));
ax=gca;
ax.XLim=[2 14];
ax.YLim=[0 1e-8];

%% Compute correlation between entrainment fidelity of STIM epochs and speed variability


%% Plot the results
figure
subplot(4,3,1)
plot(out.f,mean(out.BL_spectrum),'color','k','Linewidth', 2)

subplot(4,3,2)
plot(out.f,mean(out.STIM_spectrum),'color',[0.8 0 0],'Linewidth', 2)

subplot(4,3,4)
[B,sorting_BL_idx] = sort(out.BL_peak_frequency, 'descend') % sort to dominant freq
imagesc(out.f,1:size(out.BL_spectrum,1),out.BL_spectrum(sorting_BL_idx,:))

subplot(4,3,5)
if stim_frequency ~= 0
[B,sorting_STIM_idx] = sort(out.STIM_entrainment_fidelity, 'descend') % sort to dominant freq
else
[B,sorting_STIM_idx] = sort(out.STIM_peak_frequency, 'descend') % sort to dominant freq
end    
imagesc(out.f,1:size(out.STIM_spectrum,1),out.STIM_spectrum(sorting_STIM_idx,:))

if stim_frequency ~= 0
subplot(4,3,6)
plot(out.BL_entrainment_fidelity(sorting_BL_idx),1:length(out.BL_entrainment_fidelity),'color', 'k')
hold on
plot(out.STIM_entrainment_fidelity(sorting_STIM_idx),1:length(out.STIM_entrainment_fidelity),'color','r')
ax=gca;
ax.YDir = 'reverse';
end

figure
subplot(2,2,1)
time_axis = linspace(-stim_duration./Fs,stim_duration./Fs,size(out.speed_STA,1));
plot(time_axis,out.speed_STA,'color', [0.5 0.5 0.5])
hold on
plot(time_axis,out.mean_speed_STA,'color',[0.1 0.1 0.8], 'Linewidth', 2)
plot(time_axis,out.mean_speed_STA+out.SEM_speed_STA,'--','color',[0.1 0.1 0.8])
plot(time_axis,out.mean_speed_STA-out.SEM_speed_STA,'--','color',[0.1 0.1 0.8])
title 'Speed STA'
hold off

subplot(2,2,2)
theta_BL_idx = find(out.BL_peak_frequency > 5);
theta_STIM_idx = find(out.STIM_peak_frequency > 5);

scatter(out.BL_mean_speed(theta_BL_idx),out.BL_peak_frequency(theta_BL_idx));
hold on
scatter(out.STIM_mean_speed(theta_STIM_idx),out.STIM_peak_frequency(theta_STIM_idx));


end

