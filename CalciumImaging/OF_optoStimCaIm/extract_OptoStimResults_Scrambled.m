%% Load data here (ms, behav)
workingDir = pwd;

load([workingDir filesep 'ms.mat']);
load([workingDir filesep 'behav.mat']);

%% Define ideal threshold for opto stim
[stimulatedEpochs,OpticSigThreshold] = extractOptoStimScrambled(ms,behav)

%% Extract relationship between stim and activity
[zMI, p_values, baselineActivity, stimActivity, portBL] = optoStimMI(ms, behav, stimulatedEpochs)

%% How many cells are significantly modulated?
sigModulatedCells = find(p_values < 0.05);
numSigModCells = length(sigModulatedCells);
portSigModCells = numSigModCells/length(p_values);

portExcited = sum(portBL(sigModulatedCells)>1)/numSigModCells;
portInhibited = sum(portBL(sigModulatedCells)<1)/numSigModCells;

%% Plot STA of significant cells (raw and binarized)


%% TEMP: add trace smoothing or some sort of normalization
% max_val = max(ms.RawTraces(:));
% gain = 20;
% for cell_i = 1:length(sigModulatedCells)
%     figure
%     [STA_onset_binary, STA_offset_binary] = optoStimSTA(ms, behav, sigModulatedCells(cell_i), OpticSigThreshold);
%     [STA_onset_raw, STA_offset_raw] = optoStimSTA_raw(ms, behav, sigModulatedCells(cell_i), OpticSigThreshold);
%     
%     time_vector = (1:size(STA_onset_raw,2))/30;
%     
%     on_times = STA_onset_binary==1;
%     off_times = STA_onset_binary==0;
%     
%     for stim_i = 1:size(STA_onset_raw,1)
%        plot(time_vector,STA_onset_raw(stim_i,:)*gain+max_val*stim_i, 'color',[0.5 0.5 0.5]); hold on
%        plot(time_vector(on_times(stim_i,:)),STA_onset_raw(stim_i,on_times(stim_i,:))*gain+max_val*stim_i, 'color',[0.5 0.8 0.5]);hold on
%     end
% end

save([workingDir filesep 'optoStimResults.mat'], 'numSigModCells','p_values','portBL','portExcited','portInhibited','portSigModCells','sigModulatedCells','stimActivity','baselineActivity','zMI');