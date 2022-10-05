function msPlotTransients(ms, cell)
%MSPLOTTRANSIENTS Summary of this function goes here
%   Detailed explanation goes here

%% Parameters
dt = median(diff(ms.time))/1000; % Conversion from ms to s
Fs = 1/dt;
window_ahead = 1;
window_behind = 2;
window_vec = -window_ahead:dt:window_behind+dt;

acorr_maxlag = 5; % in s

hist_bin = 20;

time_vec = ms.time/1000;
current_trace = ms.RawTraces(:,cell);
peak_idx = ms.transients{cell}.peak_idx;

%% Spike-triggered average
for peak_i = 1:length(peak_idx);
   if peak_idx(peak_i)-round(window_ahead*Fs > 1 & peak_idx(peak_i)+round(window_behind*Fs) < length(current_trace));
   trans_matrix(:,peak_i) = current_trace(peak_idx(peak_i)-round(window_ahead*Fs):peak_idx(peak_i)+round(window_behind*Fs)); 
   end
end

avg_trace = mean(trans_matrix,2);

%% Autocorrelogram

[C, lags] = xcorr(current_trace,round(acorr_maxlag*Fs), 'coeff');

%% Plotting the results
figure;
%% Plotting trace
subplot(4,5,[1 5]);
plot(time_vec,current_trace, 'black'); hold on;
for peak_i=1:length(peak_idx);
line([time_vec(peak_idx(peak_i)) time_vec(peak_idx(peak_i))], [max(current_trace) 1.2*max(current_trace)], 'color', [0 0 1])
end
ax=gca;
ax.YLim = [-0.2 max(current_trace)*1.2];
ax.XLim = [0 max(time_vec)];
ax.FontSize = 12;
xlabel 'Time (s)'

%% Plotting triggered average
subplot(4,5,[6 7 11 12]);
for peak_i=1:length(peak_idx);
plot(window_vec,trans_matrix(:,peak_i), 'color',[0.8 0.8 0.8]); hold on
end
line(window_vec,avg_trace, 'color', [1 0.4 0.4], 'linewidth', 4);
hold off
xlabel 'Time (s)'
ylabel 'dF/F'
ax=gca;
ax.XLim = [min(window_vec) max(window_vec)];
ax.FontSize = 12;

%% Plotting frequency
subplot(4,5,8);
boxplot(1./diff(ms.transients{cell}.start_time));
ylabel 'Frequency Hz'
ax=gca;
ax.FontSize = 12;
ms.Frequency = mean(1./diff(ms.transients{cell}.start_time));


%% Plotting amplitude
subplot(4,5,9);
boxplot(ms.transients{cell}.prominence);
ylabel 'Amplitude (dF/F)'
ax=gca;
ax.FontSize = 12;

%% Plotting rise time
subplot(4,5,10);
boxplot(ms.transients{cell}.rise_time);
ylabel 'Rise time (ms)'
ax=gca;
ax.FontSize = 12;

%% Plotting decay time
subplot(4,5,13);
boxplot(ms.transients{cell}.decay_time);
ylabel 'Decay time (ms)'
ax=gca;
ax.FontSize = 12;

%% Plotting max rise slope
subplot(4,5,14);
boxplot(ms.transients{cell}.max_rise_slope);
ylabel 'Max rise slope (dF.ms-1)'
ax=gca;
ax.FontSize = 12;

%% Plotting max decay slope
subplot(4,5,15);
boxplot(ms.transients{cell}.max_decay_slope);
ylabel 'Max decay slope (dF.ms-1)'
ax=gca;
ax.FontSize = 12;

%% Plotting autocorrelation
subplot(4,5,16);
plot(lags/Fs,C,'color',[0 0 0.6]);
xlabel 'Time (s)'
ylabel 'Correlation coefficient'
ax=gca;
ax.XLim = [min(lags/Fs) max(lags/Fs)];
ax.FontSize = 12;


%% Plotting half-width
subplot(4,5,18);
boxplot(ms.transients{cell}.width);
ylabel 'Half width (ms)'
ax=gca;
ax.FontSize = 12;

%% Plotting rise area
subplot(4,5,19);
boxplot(ms.transients{cell}.rise_area);
ylabel 'Rise area (dF2)'
ax=gca;
ax.FontSize = 12;

%% Plotting decay area
subplot(4,5,20);
boxplot(ms.transients{cell}.decay_area);
ylabel 'Decay area (dF2)'
ax=gca;
ax.FontSize = 12;



end

