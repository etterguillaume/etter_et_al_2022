%% Load data here
load('/Users/guillaume/Dropbox (Williams Lab)/Guillaume/analyzed_data/calcium_imaging/PV986/OF/PV986_openfield_whitenoise5s_20190116/behav.mat')
load('/Users/guillaume/Dropbox (Williams Lab)/Guillaume/analyzed_data/calcium_imaging/PV986/OF/PV986_openfield_whitenoise5s_20190116/ms.mat')

%% Set opto signal threshold here
%[OpticSigThreshold] = optoStimExtractSigThreshold(behav)
OpticSigThreshold = 0.2;

%% Extract binary vector corresponding to stimulated epochs
%% Parameters
Fs = round(1/mode(diff(ms.time/1000)));
[bFilt,aFilt] = butter(2,  2/(30/2), 'low');

StimDuration = 5; % in s
StimDurationFrames = round(StimDuration*Fs);

%% Make sure time vectors contain unique time values
[unique_behavtime, IAbehav, ICbehav]=unique(behav.time);
[unique_mstime, IAms, ICms]=unique(ms.time);

numFrames = length(unique_mstime);

%% Re-aligning calcium and behavior times
optoSignal = behav.optosignal;
InterpolatedOptoSignal = interp1(behav.time(IAbehav),optoSignal(IAbehav),ms.time(IAms));
smoothedOptoSignal = smooth(InterpolatedOptoSignal, 'sgolay');
normalizedOptoSignal = smoothedOptoSignal-min(smoothedOptoSignal);
normalizedOptoSignal = normalizedOptoSignal./max(normalizedOptoSignal);

stimulatedEpochs = zeros(size(normalizedOptoSignal));
stimulatedEpochs(normalizedOptoSignal>OpticSigThreshold) = 1;
stimulatedEpochs = logical(stimulatedEpochs);

%% Extract binary signals
for cell_i = 1:size(ms.RawTraces,2)
    binaryData(:,cell_i) = extract_binary(ms.RawTraces(IAms,cell_i),Fs,2);
end

%% Extract cross MI on baseline epochs
[BL_zMI, BL_MI, BL_p_value] = extractCrossMI(binaryData(stimulatedEpochs,:), true);

%% Extract cross MI on stim epochs
[stim_zMI, stim_MI, stim_p_value] = extractCrossMI(binaryData(~stimulatedEpochs,:), true);

%% Normalize MI values
norm_BL_MI = BL_MI./max(BL_MI);
norm_stim_MI = stim_MI./max(stim_MI);

%norm_BL_MI(max(BL_MI) == 0,:) = [];

%% Compute some sort of metric
norm_MI_port_change = norm_stim_MI./norm_BL_MI;
MI_port_change = stim_MI./BL_MI;
MI_port_change = MI_port_change(:);
MI_port_change(isnan(MI_port_change)) = [];
sorted_MI_change = sort(MI_port_change);


for i = 1:size(BL_MI,1)
    for j = 1:size(BL_MI,2)
        if i == j
            norm_BL_MI(i,j) = NaN;
            norm_stim_MI(i,j) = NaN; 
            BL_MI(i,j) = NaN;
            stim_MI(i,j) = NaN;
        end
    end
end

figure
subplot(1,2,1)
imagesc(BL_MI)
title 'Baseline'
ax=gca;
ax.CLim = [0 0.01];
colorbar
daspect([1 1 1])
colormap Viridis

subplot(1,2,2)
imagesc(stim_MI)
title 'Stimulation'
ax=gca;
ax.CLim = [0 0.01];
colorbar
daspect([1 1 1])
colormap Viridis

norm_BL_MI = norm_BL_MI(:);
norm_stim_MI = norm_stim_MI(:);

norm_BL_MI(isnan(norm_BL_MI)) = [];
norm_stim_MI(isnan(norm_stim_MI)) = [];
