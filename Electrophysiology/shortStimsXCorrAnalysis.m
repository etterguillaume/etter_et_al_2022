function [out] = geShortStimsXCorrAnalysis(in)
%ge_wavcoupling Wavelet- and Hilbert-based coupling function. This filters and extracts the phase of theta oscillations,
%and find associated higher (modulated) frequencies. These are split into
%low, medium and high gamma. Data is expressed as standard-deviation from a uniform distribution of theta phases.
%Guillaume Etter: guillaume.etter@mail.mcgill.ca (contact if any issue/bug)

%% Get user input
channel1=input('Please enter the first channel to be analyzed: ');
channel2=input('Please enter the second channel to be analyzed: ');
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
f_n = (1/dt)/2;             % nyquist frequency
theta_band = theta_band./f_n;
[theta_B] = fir1(1000,theta_band);

%% Extracting data
display('Extracting and filtering data.....');
%data = zscore(in.LFP.Data(:,channel));
data1 = in.LFP.Data(:,channel1);
data2 = in.LFP.Data(:,channel2);

theta1 = filtfilt(theta_B,1,data1);
theta2 = filtfilt(theta_B,1,data2);

display('Hilbert transform.....');
theta_phase1 = angle(hilbert(theta1))*360/(2.0*pi);
theta_phase2 = angle(hilbert(theta2))*360/(2.0*pi);

%% Extracting event times
display('Detecting TTLs.....');
events=in.Events.TimeStamps{1};
events_int=diff(events);
events_freq=1./events_int;

stims(1)= round(events(1) * Fs);

%% Detecting stim start time
for i=1:length(events_freq)-1
    if events_freq(i) < 1 && events_freq(i+1)>=1.75 % Detection algorithm: detects sudden increases in TTL frequency ie putatively stimulated epochs
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

for i = 1:length(stims)
    BL_span = stims(i)-stim_duration:stims(i);
    stim_span = stims(i):stims(i)+stim_duration;
    full_span = stims(i)-stim_duration:stims(i)+stim_duration;
    
    %% Cross-correlation
    BL_theta1(:,i) = theta1(BL_span);
    BL_theta2(:,i) = theta2(BL_span);
    stim_theta1(:,i) = theta1(stim_span);
    stim_theta2(:,i) = theta2(stim_span);
    
    [BL_corr(:,i), lags] = xcorr(BL_theta1(:,i),BL_theta2(:,i),1/dt);
    [stim_corr(:,i), lags] = xcorr(stim_theta1(:,i),stim_theta2(:,i),1/dt);
    
    full_theta1(:,i) = theta1(full_span);
    full_theta2(:,i) = theta2(full_span);
    
    BL_theta_phase1(:,i) = theta_phase1(BL_span);
    BL_theta_phase2(:,i) = theta_phase2(BL_span);
    stim_theta_phase1(:,i) = theta_phase1(stim_span);
    stim_theta_phase2(:,i) = theta_phase2(stim_span);
    full_theta_phase1(:,i) = theta_phase1(full_span);
    full_theta_phase2(:,i) = theta_phase2(full_span);
    
        %% Phase coupling analysis
%     for deg_i = 1:length(deg_vec)
%         % Baselines
%         BL_I = find(BL_theta_phase1(:,i) >= deg_vec(deg_i) & BL_theta_phase1(:,i) <  deg_vec(deg_i)+deg_bin);
%         BL_amp_theta2(deg_i,i) = mean(BL_theta2(BL_I),'omitnan');
%         BL_I = find(BL_theta_phase2(:,i) >= deg_vec(deg_i) & BL_theta_phase2(:,i) <  deg_vec(deg_i)+deg_bin);
%         BL_amp_theta1(deg_i,i) = mean(BL_theta1(BL_I),'omitnan');
%         
%         % Stims
%         stim_I = find(stim_theta_phase1(:,i) >= deg_vec(deg_i) & stim_theta_phase1(:,i) <  deg_vec(deg_i)+deg_bin);
%         stim_amp_theta2(deg_i,i) = mean(stim_theta2(stim_I),'omitnan');
%         stim_I = find(stim_theta_phase2(:,i) >= deg_vec(deg_i) & stim_theta_phase2(:,i) <  deg_vec(deg_i)+deg_bin);
%         stim_amp_theta1(deg_i,i) = mean(stim_theta1(stim_I),'omitnan');
%     end
%     
%     % MI calculation
%     MI_BL_theta1(i)=(log(nbin)-(-sum((BL_amp_theta1(:,i)/sum(BL_amp_theta1(:,i))).*log((BL_amp_theta1(:,i)+eps/sum(BL_amp_theta1(:,i)))))))/log(nbin);
%     MI_BL_theta2(i)=(log(nbin)-(-sum((BL_amp_theta2(:,i)/sum(BL_amp_theta2(:,i))).*log((BL_amp_theta2(:,i)+eps/sum(BL_amp_theta2(:,i)))))))/log(nbin);
%     MI_stim_theta1(i)=(log(nbin)-(-sum((stim_amp_theta1(:,i)/sum(stim_amp_theta1(:,i))).*log((stim_amp_theta1(:,i)+eps/sum(stim_amp_theta1(:,i)))))))/log(nbin);
%     MI_stim_theta2(i)=(log(nbin)-(-sum((stim_amp_theta2(:,i)/sum(stim_amp_theta2(:,i))).*log((stim_amp_theta2(:,i)+eps/sum(stim_amp_theta2(:,i)))))))/log(nbin);
%     
    
end

figure;imagesc(full_theta_phase2'); colormap('HSV')

figure;imagesc(BL_corr); colormap('Viridis')
figure;imagesc(stim_corr); colormap('Viridis')



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