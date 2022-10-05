%% Load data here
load('/Users/guillaume/Dropbox (Williams Lab)/Guillaume/analyzed_data/Electrophysiology/Chrimson/FallWinter_2018/OF_8Hz_1mW_5sONOFF/PV950/data_PV950_OF_8Hz_1mW_5sONOFF.mat')

% PV950: max ripple power in ch2
% PV951: broken ground?
% PV952: max ripple power in ch5
% PV953: max ripple power in ch5
% PV1036: max ripple power in ch7
% PV1040: max ripple power in ch6
% PV1045: max ripple power in ch7

Fs = in.LFP.SamplingRate;
data = in.LFP.Data(:,7);
events = in.Events.TimeStamps{1};

%% 1 - Identify theta/delta ratio to only consider periods of inactivity
[theta_B] = fir1(1000,[4 12]/(Fs/2));
[delta_B] = fir1(1000,[1 4]/(Fs/2));
theta = filtfilt(theta_B,1,data);
delta = filtfilt(delta_B,1,data);
theta_power = abs(hilbert(theta));
delta_power = abs(hilbert(delta));

TDR = zscore(theta_power./delta_power);
LIA_state = TDR < 0;


%% 2 - Isolate periods of inactivity




%% 3 - Filter ripple frequency



%% 4 - Z-score



%% 5 - Count ripple peaks (peak detection?)




