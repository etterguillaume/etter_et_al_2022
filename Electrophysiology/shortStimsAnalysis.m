function [out] = geShortStimsAnalysis(in)
%ge_wavcoupling Wavelet- and Hilbert-based coupling function. This filters and extracts the phase of theta oscillations,
%and find associated higher (modulated) frequencies. These are split into
%low, medium and high gamma. Data is expressed as standard-deviation from a uniform distribution of theta phases.
%Guillaume Etter: guillaume.etter@mail.mcgill.ca (contact if any issue/bug)

%% Get user input
channel=input('Please enter the channel to be analyzed: ');
stim_duration = input('Please indicate the duration of stimulations (eg 2s): ');
stim_frequency = input('Please indicate the frequency of stimulations (eg 40): ');
start=in.Spread(1);
stop=in.Spread(2);

%% Parameters
Fs = in.LFP.SamplingRate;
dt = in.LFP.dt;

t = (0:dt:stim_duration)';     % seconds                   
stim_duration = round(stim_duration * Fs);

deg_bin = 9;
deg_vec = -180:deg_bin:180-deg_bin;
bin_center = deg_bin/2;
deg_vec_plot = deg_vec + bin_center;
nbin = length(deg_vec);

%% Chronux parameters
params.Fs = Fs;            % sampling frequency
params.tapers = [1 1];    % taper parameter TODO double check p 12!!! TW K default[1 1]
params.pad = 4;           %
params.fpass = [1 300];  % frequency range ; default [3 55]
params.err = 0;           % 0 for no error calculated
params.trialave = 0;      % when true average the trials

% Filter parameters
theta_band = [4 12];
slow_gamma_band = [30 59];
fast_gamma_band = [61 120];

f_n = (1/dt)/2;             % nyquist frequency

theta_band = theta_band./f_n;
slow_gamma_band = slow_gamma_band./f_n;
fast_gamma_band = fast_gamma_band./f_n;

[theta_B] = fir1(1000,theta_band);
[slow_gamma_B] = fir1(1000,slow_gamma_band);
[fast_gamma_B] = fir1(1000,fast_gamma_band);

%% Extracting data
display('Extracting and filtering data.....');
%data = zscore(in.LFP.Data(:,channel));
data = in.LFP.Data(:,channel);
theta = filtfilt(theta_B,1,data);
sgamma = filtfilt(slow_gamma_B,1,data);
fgamma = filtfilt(fast_gamma_B,1,data);

display('Hilbert transform.....');
theta_power = abs(hilbert(theta));
%theta_phase = ge_get_theta_phase(data,Fs);
%theta_phase = theta_phase - 180;
theta_phase = angle(hilbert(theta))*360/(2.0*pi);
sgamma_power = abs(hilbert(sgamma));
fgamma_power = abs(hilbert(fgamma));


%% Extracting event times
display('Detecting TTLs.....');
events=in.Events.TimeStamps{1};
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
if stims(end)+ stim_duration >= round(stop * Fs)
    stims(end)=[];
end

if stims(1)- stim_duration <= round(start * Fs)
    stims(1)=[];
end

%% Extracting data segments
display('Extracting data segments.....');

for i = 1:length(stims);
    BL_span = stims(i)-stim_duration:stims(i);
    stim_span = stims(i):stims(i)+stim_duration;
    
    [BL_spectrum(:,i),f] = mtspectrumc(data(BL_span),params);
    [stim_spectrum(:,i),f] = mtspectrumc(data(stim_span),params);
    
    % Normalizing control signals
    data_floor = min(data(BL_span));
    data_ceil = max(data(BL_span));
    control_sig = (data_ceil-data_floor).*sin(2*pi*stim_frequency*t)+data_floor; % Generate control sine wave
    control_noise = (data_ceil-data_floor).*rand(stim_duration+1,1)+data_floor; % Generate white noise
    
    [control_sig_spectrum(:,i),f] = mtspectrumc(control_sig,params);
    [control_noise_spectrum(:,i),f] = mtspectrumc(control_noise,params);
    
    %% OS strength analysis
    % Normalize spectra
    BL_norm_theta_spectr(:,i) = BL_spectrum(:,i)-min(BL_spectrum(f > 4 & f < 12,i));
    stim_norm_theta_spectr(:,i) = stim_spectrum(:,i)-min(stim_spectrum(f > 4 & f < 12,i));
    control_sig_norm_theta_spectr(:,i) = control_sig_spectrum(:,i) - min(control_sig_spectrum(f > 4 & f < 12,i));
    control_noise_norm_theta_spectr(:,i) = control_noise_spectrum(:,i) - min(control_noise_spectrum(f > 4 & f < 12,i));
    
    out.BL_theta_OS(:,i) = sum(BL_norm_theta_spectr(f > stim_frequency - 2 & f < stim_frequency + 2,i))/sum(BL_norm_theta_spectr(f > 4 & f < 12,i));
    out.stim_theta_OS(:,i) = sum(stim_norm_theta_spectr(f > stim_frequency - 2 & f < stim_frequency + 2,i))/sum(stim_norm_theta_spectr(f > 4 & f < 12,i));
    out.control_sig_theta_OS(:,i) = sum(control_sig_norm_theta_spectr(f > stim_frequency - 2 & f < stim_frequency + 2,i))/sum(control_sig_norm_theta_spectr(f > 4 & f < 12,i));
    out.control_noise_theta_OS(:,i) = sum(control_noise_norm_theta_spectr(f > stim_frequency - 2 & f < stim_frequency + 2,i))/sum(control_noise_norm_theta_spectr(f > 4 & f < 12,i));
    
    %% Peak/band powers
    BL_theta(:,i) = theta(BL_span);
    BL_theta_phase(:,i) = theta_phase(BL_span);
    BL_sgamma_power(:,i) = sgamma_power(BL_span);
    BL_fgamma_power(:,i) = fgamma_power(BL_span);
    
    out.BL_stimfreq_power(:,i) = mean(BL_spectrum(f > stim_frequency - 0.5 & f < stim_frequency + 0.5,i));
    out.stim_stimfreq_power(:,i) = mean(stim_spectrum(f > stim_frequency - 0.5 & f < stim_frequency + 0.5,i));
    
    out.BL_min_theta_power(:,i) = min(BL_spectrum(f>4 & f<12));
    out.stim_min_theta_power(:,i) = min(stim_spectrum(f>4 & f<12));
    
    out.control_sig_min_theta_power(:,i) = min(control_sig_spectrum(f>4 & f<12,i));
    out.control_noise_min_theta_power(:,i) = min(control_noise_spectrum(f>4 & f<12,i));
    
    %% Band power
    out.BL_delta_FFT_power(:,i) = mean(BL_spectrum(f>1 & f<4,i));
    out.BL_theta_FFT_power(:,i) = mean(BL_spectrum(f>4 & f<12,i));
    out.BL_alphabeta_FFT_power(:,i) = mean(BL_spectrum(f>12 & f<30,i));
    out.BL_sgamma_FFT_power(:,i) = mean(BL_spectrum(f>30 & f<59,i));
    out.BL_fgamma_FFT_power(:,i) = mean(BL_spectrum(f>61 & f<120,i));
    out.BL_HFO_FFT_power(:,i) = mean(BL_spectrum(f>120 & f<250,i));
    
    out.stims_delta_FFT_power(:,i) = mean(stim_spectrum(f>1 & f<4,i));
    out.stims_theta_FFT_power(:,i) = mean(stim_spectrum(f>4 & f<12,i));
    out.stims_alphabeta_FFT_power(:,i) = mean(stim_spectrum(f>12 & f<30,i));
    out.stims_sgamma_FFT_power(:,i) = mean(stim_spectrum(f>30 & f<59,i));
    out.stims_fgamma_FFT_power(:,i) = mean(stim_spectrum(f>61 & f<120,i));
    out.stims_HFO_FFT_power(:,i) = mean(stim_spectrum(f>120 & f<250,i));
    
    out.control_sig_theta_FFT_power(:,i) = mean(control_sig_spectrum(f>4 & f<12,i));
    out.control_noise_theta_FFT_power(:,i) = mean(control_noise_spectrum(f>4 & f<12,i));
    
    %% Peak power
    out.BL_delta_peak_power(:,i) = max(BL_spectrum(f>1 & f<4,i));
    out.BL_theta_peak_power(:,i) = max(BL_spectrum(f>4 & f<12,i));
    out.BL_alphabeta_peak_power(:,i) = max(BL_spectrum(f>12 & f<30,i));
    out.BL_sgamma_peak_power(:,i) = max(BL_spectrum(f>30 & f<59,i));
    out.BL_fgamma_peak_power(:,i) = max(BL_spectrum(f>61 & f<120,i));
    out.BL_HFO_peak_power(:,i) = max(BL_spectrum(f>120 & f<250,i));
    
    out.stims_delta_peak_power(:,i) = max(stim_spectrum(f>1 & f<4,i));
    out.stims_theta_peak_power(:,i) = max(stim_spectrum(f>4 & f<12,i));
    out.stims_alphabeta_peak_power(:,i) = max(stim_spectrum(f>12 & f<30,i));
    out.stims_sgamma_peak_power(:,i) = max(stim_spectrum(f>30 & f<59,i));
    out.stims_fgamma_peak_power(:,i) = max(stim_spectrum(f>61 & f<120,i));
    out.stims_HFO_peak_power(:,i) = max(stim_spectrum(f>120 & f<250,i));
    
    out.control_sig_peak_theta_power(:,i) = max(control_sig_spectrum(f>4 & f<12,i));
    out.control_noise_peak_theta_power(:,i) = max(control_noise_spectrum(f>4 & f<12,i));
    
    %% Hilbert analysis
    stims_theta(:,i) = theta(stim_span);
    stims_theta_phase(:,i) = theta_phase(stim_span);
    stims_sgamma_power(:,i) = sgamma_power(stim_span);
    stims_fgamma_power(:,i) = fgamma_power(stim_span);
    
    BL_theta_power_peak(i) = max(theta_power(BL_span));
    BL_sgamma_power_peak(i) = max(sgamma_power(BL_span));
    BL_fgamma_power_peak(i) = max(fgamma_power(BL_span));
    
    stims_theta_power_peak(i) = max(theta_power(stim_span));
    stims_sgamma_power_peak(i) = max(sgamma_power(stim_span));
    stims_fgamma_power_peak(i) = max(fgamma_power(stim_span));
    
    %% Phase coupling analysis
    for deg_i = 1:length(deg_vec)
        % Baselines
        BL_I = find(BL_theta_phase(:,i) >= deg_vec(deg_i) & BL_theta_phase(:,i) <  deg_vec(deg_i)+deg_bin);
        BL_amp_theta(deg_i,i) = mean(BL_theta(BL_I),'omitnan');
        BL_Meansgamma_Amp(deg_i,i)=mean(BL_sgamma_power(BL_I,i),'omitnan');
        BL_Meanfgamma_Amp(deg_i,i)=mean(BL_fgamma_power(BL_I,i),'omitnan');
        
        % Stims
        stims_I = find(stims_theta_phase(:,i) >= deg_vec(deg_i) & stims_theta_phase(:,i) <  deg_vec(deg_i)+deg_bin);
        stims_amp_theta(deg_i,i) = mean(stims_theta(stims_I),'omitnan');
        stims_Meansgamma_Amp(deg_i,i)=mean(stims_sgamma_power(stims_I,i),'omitnan');
        stims_Meanfgamma_Amp(deg_i,i)=mean(stims_fgamma_power(stims_I,i),'omitnan');
    end
    
    % MI calculation
    MI_BL_sgamma(i)=(log(nbin)-(-sum((BL_Meansgamma_Amp(:,i)/sum(BL_Meansgamma_Amp(:,i))).*log((BL_Meansgamma_Amp(:,i)+eps/sum(BL_Meansgamma_Amp(:,i)))))))/log(nbin);
    MI_BL_fgamma(i)=(log(nbin)-(-sum((BL_Meanfgamma_Amp(:,i)/sum(BL_Meanfgamma_Amp(:,i))).*log((BL_Meanfgamma_Amp(:,i)+eps/sum(BL_Meanfgamma_Amp(:,i)))))))/log(nbin);
    MI_stims_sgamma(i)=(log(nbin)-(-sum((stims_Meansgamma_Amp(:,i)/sum(stims_Meansgamma_Amp(:,i))).*log((stims_Meansgamma_Amp(:,i)+eps/sum(stims_Meansgamma_Amp(:,i)))))))/log(nbin);
    MI_stims_fgamma(i)=(log(nbin)-(-sum((stims_Meanfgamma_Amp(:,i)/sum(stims_Meanfgamma_Amp(:,i))).*log((stims_Meanfgamma_Amp(:,i)+eps/sum(stims_Meanfgamma_Amp(:,i)))))))/log(nbin);

end

    %% Oscillation strength
%     out.BL_OS=out.BL_theta_peak_power-out.BL_min_theta_power;
%     out.stim_OS=out.stims_theta_peak_power-out.stim_min_theta_power;
%     out.control_sig_OS=out.control_sig_peak_theta_power-out.control_sig_min_theta_power;
%     out.control_noise_OS=out.control_noise_peak_theta_power-out.control_noise_min_theta_power;
%     
%     out.BL_theta_OS=(out.BL_OS-out.control_noise_OS)./out.control_sig_OS;
%     out.stim_theta_OS=(out.stim_OS-out.control_noise_OS)./out.control_sig_OS;
    
    %% Normalize values here
    %% Band power
    BL_theta_band = mean(out.BL_theta_FFT_power);
    out.BL_delta_FFT_power = mean(out.BL_delta_FFT_power)/BL_theta_band;
    out.BL_theta_FFT_power = mean(out.BL_theta_FFT_power)/BL_theta_band;
    out.BL_alphabeta_FFT_power = mean(out.BL_alphabeta_FFT_power)/BL_theta_band;
    out.BL_sgamma_FFT_power = mean(out.BL_sgamma_FFT_power)/BL_theta_band;
    out.BL_fgamma_FFT_power = mean(out.BL_fgamma_FFT_power)/BL_theta_band;
    out.BL_HFO_FFT_power = mean(out.BL_HFO_FFT_power)/BL_theta_band;
    out.stims_delta_FFT_power = mean(out.stims_delta_FFT_power)/BL_theta_band;
    out.stims_theta_FFT_power = mean(out.stims_theta_FFT_power)/BL_theta_band;
    out.stims_alphabeta_FFT_power = mean(out.stims_alphabeta_FFT_power)/BL_theta_band;
    out.stims_sgamma_FFT_power = mean(out.stims_sgamma_FFT_power)/BL_theta_band;
    out.stims_fgamma_FFT_power = mean(out.stims_fgamma_FFT_power)/BL_theta_band;
    out.stims_HFO_FFT_power = mean(out.stims_HFO_FFT_power)/BL_theta_band;
    
    %% Peak power
    BL_theta_peak = mean(out.BL_stimfreq_power);
    %BL_theta_peak = mean(out.BL_theta_peak_power);
    out.BL_delta_peak_power = mean(out.BL_delta_peak_power)/BL_theta_peak;
    %out.BL_theta_peak_power = mean(out.BL_theta_peak_power)/BL_theta_peak;
    out.BL_theta_peak_power = mean(out.BL_stimfreq_power)/BL_theta_peak;
    out.BL_alphabeta_peak_power = mean(out.BL_alphabeta_peak_power)/BL_theta_peak;
    out.BL_sgamma_peak_power = mean(out.BL_sgamma_peak_power)/BL_theta_peak;
    out.BL_fgamma_peak_power = mean(out.BL_fgamma_peak_power)/BL_theta_peak;
    out.BL_HFO_peak_power = mean(out.BL_HFO_peak_power)/BL_theta_peak;
    out.stims_delta_peak_power = mean(out.stims_delta_peak_power)/BL_theta_peak;
    %out.stims_theta_peak_power = mean(out.stims_theta_peak_power)/BL_theta_peak;
    out.stims_theta_peak_power = mean(out.stim_stimfreq_power)/BL_theta_peak;
    out.stims_alphabeta_peak_power = mean(out.stims_alphabeta_peak_power)/BL_theta_peak;
    out.stims_sgamma_peak_power = mean(out.stims_sgamma_peak_power)/BL_theta_peak;
    out.stims_fgamma_peak_power = mean(out.stims_fgamma_peak_power)/BL_theta_peak;
    out.stims_HFO_peak_power = mean(out.stims_HFO_peak_power)/BL_theta_peak;
    
%% Ratios    
theta_power_ratio = mean(stims_theta_power_peak)./mean(BL_theta_power_peak);
sgamma_power_ratio = mean(stims_sgamma_power_peak)./mean(BL_sgamma_power_peak);
fgamma_power_ratio = mean(stims_fgamma_power_peak)./mean(BL_fgamma_power_peak);

MI_ratio_sgamma = mean(MI_stims_sgamma)./mean(MI_BL_sgamma);
MI_ratio_fgamma = mean(MI_stims_fgamma)./mean(MI_BL_fgamma);

%% Averageing for display
S_BL_mean = mean(BL_spectrum,2);
S_BL_SEM = std(BL_spectrum,[],2)./sqrt(size(BL_spectrum,2));
S_stim_mean = mean(stim_spectrum,2);
S_stim_SEM = std(stim_spectrum,[],2)./sqrt(size(stim_spectrum,2));


BL_theta_amp_mean = mean(BL_amp_theta,2,'omitnan');
BL_theta_amp_SEM = std(BL_amp_theta,0,2,'omitnan')./sqrt(length(stims));

BL_sgamma_PAC_mean = mean(BL_Meansgamma_Amp,2,'omitnan');
BL_sgamma_PAC_SEM = std(BL_Meansgamma_Amp,0,2,'omitnan')./sqrt(length(stims));

BL_fgamma_PAC_mean = mean(BL_Meanfgamma_Amp,2,'omitnan');
BL_fgamma_PAC_SEM = std(BL_Meanfgamma_Amp,0,2,'omitnan')./sqrt(length(stims));

stim_theta_amp_mean = mean(stims_amp_theta,2,'omitnan');
stim_theta_amp_SEM = std(stims_amp_theta,0,2,'omitnan')./sqrt(length(stims));

stim_sgamma_PAC_mean = mean(stims_Meansgamma_Amp,2,'omitnan');
stim_sgamma_PAC_SEM = std(stims_Meansgamma_Amp,0,2,'omitnan')./sqrt(length(stims));

stim_fgamma_PAC_mean = mean(stims_Meanfgamma_Amp,2,'omitnan');
stim_fgamma_PAC_SEM = std(stims_Meanfgamma_Amp,0,2,'omitnan')./sqrt(length(stims));

theta_phase_ratio = stim_theta_amp_mean./BL_theta_amp_mean;
sgamma_phase_ratio = stim_sgamma_PAC_mean./BL_sgamma_PAC_mean;
fgamma_phase_ratio = stim_fgamma_PAC_mean./BL_fgamma_PAC_mean;

%% Sorting spectra
[~, max_BL_spectra_idx] = max(BL_spectrum);
[~, sorted_BL_spectra_idx] = sort(max_BL_spectra_idx,'descend');
[max_stim_spectra_idx] = max(stim_spectrum);
[~, sorted_stim_spectra_idx] = sort(max_stim_spectra_idx,'descend');

%% Normalizing spectra
max_baseline_spectr = max(BL_spectrum(:));
max_baseline_mean = max(mean(BL_spectrum,2));
norm_BL_spectrum = BL_spectrum./max_baseline_spectr*100;
norm_stim_spectrum = stim_spectrum./max_baseline_spectr*100;
norm_BL_mean = mean(BL_spectrum,2)./max_baseline_mean*100;
norm_stim_mean = mean(stim_spectrum,2)./max_baseline_mean*100;

%% Plotting the results
figure
subplot(2,2,1)
plot(f,mean(norm_BL_mean,2),'color','k')
ax1 = gca;
ax1.XLim = [2 14];
ax1.YLim = [0 100];

subplot(2,2,2)
plot(f,mean(norm_stim_mean,2),'color', [0.8 0 0])
ax2=gca;
ax2.XLim = [2 14];
ax2.YLim = [0 4500];

subplot(2,2,3)
imagesc(f,1:length(stims),norm_BL_spectrum(:,sorted_BL_spectra_idx)');
colormap viridis
ax3 = gca;
ax3.XLim = [2 14];
ax3.CLim = [0 100];
colorbar

subplot(2,2,4)
imagesc(f,1:length(stims),norm_stim_spectrum(:,sorted_stim_spectra_idx)');
colormap viridis
ax4 = gca;
ax4.XLim = [2 14];
ax4.CLim = [0 800];
colorbar

linkaxes([ax1 ax2 ax3 ax4], 'x')

figure
plot(f,S_BL_mean, 'color', 'black', 'LineWidth', 2); hold on;
plot(f,S_BL_mean+S_BL_SEM, 'color', 'black', 'LineStyle', '--');
plot(f,S_BL_mean-S_BL_SEM, 'color', 'black', 'LineStyle', '--');
plot(f,S_stim_mean, 'color', 'red', 'LineWidth', 2); hold on;
plot(f,S_stim_mean+S_stim_SEM, 'color', 'red', 'LineStyle', '--');
plot(f,S_stim_mean-S_stim_SEM, 'color', 'red', 'LineStyle', '--');
line([4 4], ylim, 'LineStyle', '--');
line([12 12], ylim, 'LineStyle', '--');
ax=gca;
ax.XLim = [1 120];
ax.XScale = 'log';
ax.YScale = 'linear';
ax.XTick = [1 2 3 4 5 6 7 8 9 10 12 14 16 18 20 25 30 40 50 60 70 80 90 100];
ax.FontSize = 18;
title 'Spectrum';
xlabel 'Frequency (Hz)';
ylabel 'Power (mV2.Hz-1)'; hold off;

figure
subplot(3,1,1);
plot(deg_vec_plot, BL_theta_amp_mean, 'color', 'black','LineWidth',2); hold on
plot(deg_vec_plot, BL_theta_amp_mean+BL_theta_amp_SEM, 'color', 'black');
plot(deg_vec_plot, BL_theta_amp_mean-BL_theta_amp_SEM, 'color', 'black');
plot(deg_vec_plot, stim_theta_amp_mean, 'color', 'red','LineWidth',2);
plot(deg_vec_plot, stim_theta_amp_mean+stim_theta_amp_SEM, 'color', 'red');
plot(deg_vec_plot, stim_theta_amp_mean-stim_theta_amp_SEM, 'color', 'red');
hold off
ylabel 'Theta amplitude (SD)'
ax=gca;
ax.FontSize = 18;
ax.XTick = [-180 -90 0 90 180];

subplot(3,1,2);
plot(deg_vec_plot, BL_sgamma_PAC_mean, 'color', 'black','LineWidth',2); hold on
plot(deg_vec_plot, BL_sgamma_PAC_mean+BL_sgamma_PAC_SEM, 'color', 'black');
plot(deg_vec_plot, BL_sgamma_PAC_mean-BL_sgamma_PAC_SEM, 'color', 'black');

plot(deg_vec_plot, stim_sgamma_PAC_mean, 'color', 'red','LineWidth',2);
plot(deg_vec_plot, stim_sgamma_PAC_mean+stim_sgamma_PAC_SEM, 'color', 'red');
plot(deg_vec_plot, stim_sgamma_PAC_mean-stim_sgamma_PAC_SEM, 'color', 'red');

hold off
ylabel 'Slow gamma amplitude (SD)'
ax=gca;
ax.FontSize = 18;
ax.XTick = [-180 -90 0 90 180];

subplot(3,1,3);
plot(deg_vec_plot, BL_fgamma_PAC_mean, 'color', 'black','LineWidth',2); hold on
plot(deg_vec_plot, BL_fgamma_PAC_mean+BL_fgamma_PAC_SEM, 'color', 'black'); hold on
plot(deg_vec_plot, BL_fgamma_PAC_mean-BL_fgamma_PAC_SEM, 'color', 'black'); hold on

plot(deg_vec_plot, stim_fgamma_PAC_mean, 'color', 'red','LineWidth',2);
plot(deg_vec_plot, stim_fgamma_PAC_mean+stim_fgamma_PAC_SEM, 'color', 'red');
plot(deg_vec_plot, stim_fgamma_PAC_mean-stim_fgamma_PAC_SEM, 'color', 'red');
hold off
ylabel 'Fast gamma amplitude (SD)'
xlabel 'Theta phase (º)'
ax=gca;
ax.FontSize = 18;
ax.XTick = [-180 -90 0 90 180];

%% Formating the output
out.Location = in.Location;
out.Channel = in.LFP.ChannelNames(channel);
out.Spread = in.Spread;
out.SamplingRate = in.LFP.SamplingRate;
out.dt = in.LFP.dt;
out.Phases = deg_vec;

    out.BL_spectrum = BL_spectrum;
    out.stim_spectrum = stim_spectrum;

    out.StimFreq.MeanBLPower = mean(sum(BL_spectrum(f < stim_frequency +0.5 & f > stim_frequency - 0.5),2));
    out.StimFreq.MeanStimPower = mean(sum(stim_spectrum(f < stim_frequency +0.5 & f > stim_frequency - 0.5),2));

    out.Theta.amp_BL = BL_theta_amp_mean;
    out.Theta.amp_BL_SEM = BL_theta_amp_SEM;
    out.Theta.amp_stim = stim_theta_amp_mean;
    out.Theta.amp_stim_SEM = stim_theta_amp_SEM;
    out.Theta.PeakPower_BL = BL_theta_power_peak;
    out.Theta.PeakPower_stims = stims_theta_power_peak;
    out.Theta.PeakPower_ratio = theta_power_ratio;
    out.Theta.PhasePower = theta_phase_ratio;
    
    out.Sgamma.BL = BL_sgamma_PAC_mean;
    out.Sgamma.BL_SEM = BL_sgamma_PAC_SEM;
    out.Sgamma.stim = stim_sgamma_PAC_mean;
    out.Sgamma.stim_SEM = stim_sgamma_PAC_SEM;
    out.Sgamma.PeakPower_BL = BL_sgamma_power_peak;
    out.Sgamma.PeakPower_stims = stims_sgamma_power_peak;
    out.Sgamma.PeakPower_ratio = sgamma_power_ratio;
    out.Sgamma.MI_BL = MI_BL_sgamma;
    out.Sgamma.MI_stims = MI_stims_sgamma;
    out.Sgamma.MI_ratio = MI_ratio_sgamma;
    out.Sgamma.PhasePower = sgamma_phase_ratio;
    
    out.Fgamma.BL = BL_fgamma_PAC_mean;
    out.Fgamma.BL_SEM = BL_fgamma_PAC_SEM;
    out.Fgamma.stim = stim_fgamma_PAC_mean;
    out.Fgamma.stim_SEM = stim_fgamma_PAC_SEM;
    out.Fgamma.PeakPower_BL = BL_fgamma_power_peak;
    out.Fgamma.PeakPower_stims = stims_fgamma_power_peak;
    out.Fgamma.PeakPower_ratio = fgamma_power_ratio;
    out.Fgamma.MI_BL = MI_BL_fgamma;
    out.Fgamma.MI_stims = MI_stims_fgamma;
    out.Fgamma.MI_ratio = MI_ratio_fgamma;
    out.Fgamma.PhasePower = fgamma_phase_ratio;
    
    out.N_stims = length(stims);

end