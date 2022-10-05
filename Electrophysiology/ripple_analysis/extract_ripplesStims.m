function [ out ] = extract_ripplesStims( in )
%GE_RIPPLES Summary of this function goes here
%   Detailed explanation goes here
warning off all

channel=input('Please enter the channel to be analyzed: ');

stim_duration = input('Please indicate the duration of stimulations (eg 2s): ');

start=in.Spread(1);
stop=in.Spread(2);

%% Parameters
Fs = in.LFP.SamplingRate;
dt = in.LFP.dt;

stim_duration = round(stim_duration * Fs);

%% Filter parameters
 f_band = [150 250];         % Ripple frequency band
 N = 1000;                   % order
 f_n = (1/dt)/2;             % nyquist frequency
 f_band = f_band./f_n;       % convert band pass to 0 to 1 value where 1 = f_n
 
[B] = fir1(N,f_band);

MinPeakHeight = 4;

minpeakwidth = 0;
minpeakdistance = 0.03;

%% Extracting event times
display('Detecting TTLs.....');
events=in.Events.TimeStamps{1};
events_int=diff(events);
events_freq=1./events_int;

stims(1)= round(events(1) * Fs);

%% Detecting stim start time
for i=1:length(events_freq)-1
    if events_freq(i) < 1 && events_freq(i+1)>=1.75 % Detection algorithm: detects sudden increases in TTL frequency ie putatively stimulated epochs
        stims(end+1) =round(events(i+1) * Fs);
    end
end

%% Removing stims that overlap with recording segment edges
if stims(end)+ stim_duration >= round(stop * Fs)
    stims(end)=[];
end

if stims(1)- stim_duration <= round(start * Fs)
    stims(1)=[];
end

%% Extracting input data from the structure
data = in.LFP.Data(:,channel);

ripple_sig = filtfilt(B,1,data); % Filter the data in the range specified by f_bands
ripple_power = abs(hilbert(ripple_sig));
ripple_Z = zscore(ripple_power);
    
    avg_BL_peak_pwr =[];
    num_BL_ripples = [];
    avg_BL_width = [];
    avg_BL_prominence = [];
    avg_stim_peak_pwr = [];
    num_stim_ripples = [];
    avg_stim_width = [];
    avg_stim_prominence = [];

for i = 1:length(stims);    
    %% Find segments
    BL_span = stims(i)-stim_duration:stims(i);
    stim_span = stims(i):stims(i)+stim_duration;
    
    [BL_peak_pwr,BL_peak_loc, BL_width, BL_prominence]=findpeaks(ripple_Z(BL_span),Fs,'MinPeakHeight',MinPeakHeight, 'MinPeakWidth', minpeakwidth, 'MinPeakDistance',minpeakdistance);
    [stim_peak_pwr,stim_peak_loc, stim_width, stim_prominence]=findpeaks(ripple_Z(stim_span),Fs,'MinPeakHeight',MinPeakHeight, 'MinPeakWidth', minpeakwidth, 'MinPeakDistance',minpeakdistance);

    if ~isempty(BL_peak_loc)
    avg_BL_peak_pwr(end+1) = mean(BL_peak_pwr);
    num_BL_ripples(end+1) = numel(BL_peak_loc);
    avg_BL_width(end+1) = mean(BL_width);
    avg_BL_prominence(end+1) = mean(BL_prominence);
    else
    avg_BL_peak_pwr(end+1) = NaN;
    num_BL_ripples(end+1) = 0;
    avg_BL_width(end+1) = NaN;
    avg_BL_prominence(end+1) = NaN;
    end
    
    if ~isempty(stim_peak_loc)
    avg_stim_peak_pwr(end+1) = mean(stim_peak_pwr);
    num_stim_ripples(end+1) = numel(stim_peak_loc);
    avg_stim_width(end+1) = mean(stim_width);
    avg_stim_prominence(end+1) = mean(stim_prominence);
    else
    avg_stim_peak_pwr(end+1) = NaN;
    num_stim_ripples(end+1) = 0;
    avg_stim_width(end+1) = NaN;
    avg_stim_prominence(end+1) = NaN;
    end
end

%% Analysis
out.numStims = length(stims);

out.BL_num_ripples_mean = mean(num_BL_ripples,'omitnan');
out.BL_num_ripples_SEM = std(num_BL_ripples,[],'omitnan')./sqrt(out.numStims);
out.BL_ripple_power_mean = mean(avg_BL_peak_pwr,'omitnan');
out.BL_ripple_power_SEM = std(avg_BL_peak_pwr,[],'omitnan')./sqrt(out.numStims);
out.BL_ripple_width_mean = mean(avg_BL_width,'omitnan');
out.BL_ripple_width_SEM = std(avg_BL_width,[],'omitnan')./sqrt(out.numStims);

out.stim_num_ripples_mean = mean(num_stim_ripples,'omitnan');
out.stim_num_ripples_SEM = std(num_stim_ripples,[],'omitnan')./sqrt(out.numStims);
out.stim_ripple_power_mean = mean(avg_stim_peak_pwr,'omitnan');
out.stim_ripple_power_SEM = std(avg_stim_peak_pwr,[],'omitnan')./sqrt(out.numStims);
out.stim_ripple_width_mean = mean(avg_stim_width,'omitnan');
out.stim_ripple_width_SEM = std(avg_stim_width,[],'omitnan')./sqrt(out.numStims);

end

