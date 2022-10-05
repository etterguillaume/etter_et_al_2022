function [ out ] = plot_ripples( in )
%GE_RIPPLES Summary of this function goes here
%   Detailed explanation goes here

channel=input('Please enter the channel to be analyzed: ');

disp('Length of the current experiment: ');
disp(in.Spread(2));

start = input('Please indicate the starting time of the spectrogram (in s): ');    
if start == 0;
   start=in.Spread(1);
end

if isempty(start);
   start=in.Spread(1);
end

stop = input('Please indicate the ending time of the spectrogram (in s). Leave empty if you want to display until the end of the recording: ');   
if stop == 0;
   stop=in.Spread(2);
end
if isempty(stop);
   stop=in.Spread(2);
end

%% Parameters
Fs = in.LFP.SamplingRate;
dt = in.LFP.dt;
dt_start = start * Fs;
dt_stop = stop * Fs;

%% Filter parameters
 f_band = [150 250];
 N = 1000;                   % order
 f_n = (1/dt)/2;             % nyquist frequency
 f_band = f_band./f_n;       % convert band pass to 0 to 1 value where 1 = f_n
 
[B] = fir1(N,f_band);

SD_threshold = 4;
minpeakwidth = 0;
minpeakdistance = 0.03;

%% Extracting input data from the structure
span = dt_start:dt_stop;
SPAN = start:dt:stop;

for i= 1:length(span);
    data(i,:) = in.LFP.Data(span(i),channel);
end

ripple_sig = filtfilt(B,1,data); % Filter the data in the range specified by f_bands

ripple_power = abs(hilbert(ripple_sig));

ripple_Z = zscore(ripple_power);


%% Ripple peak detection
[peak_pwr,peak_loc, width, prominence]=findpeaks(ripple_Z,Fs,'MinPeakHeight',SD_threshold, 'MinPeakWidth', minpeakwidth, 'MinPeakDistance',minpeakdistance);

%% Ripples analysis

ripples_ISI=diff(peak_loc);

ripples_inst_freq = 1./ripples_ISI;

%% plot the results:

if SPAN(end)-SPAN(1)<1800;
figure
subplot(3,1,1);
plot(SPAN, data, 'color', 'black');
ax1=gca;

subplot(3,1,2);
plot(SPAN,ripple_sig)
ax2=gca;

subplot(3,1,3);
plot(SPAN, ripple_Z, 'color', 'black'); hold on;
line(xlim, [SD_threshold SD_threshold], 'color','r','LineWidth', 4);

hold on;
plot(SPAN(round(peak_loc*Fs)),peak_pwr, 'v', 'color', 'r');
ax3=gca;
hold off

linkaxes([ax1 ax2 ax3],'x');

end

%% Formatting the output

out.Ripples_loc = peak_loc;
out.Ripples_width = width;
out.Ripples_avg_width = mean(width);
out.Ripples_prominence = prominence;
out.Ripples_ISI = ripples_ISI;
out.Ripples_inst_freq = ripples_inst_freq;
%out.Ripples_epochs = ripple_epochs;

end

