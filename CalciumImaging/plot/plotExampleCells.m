function plot_CellTypes
%PLOT_CELLTYPES Summary of this function goes here
%   Detailed explanation goes here

mice = {'PV986', 'PV988', 'PV989', 'PV990', 'PV991'};
session2use = {'LMT1_20190312', 'LMT2_20190313', 'LMT5_20190318', 'LMT6_20190319'};

pooled_spatial_pvalue = [];
pooled_spatial_tuning_curves = [];
pooled_spatial_MI = [];

pooled_time_pvalue = [];
pooled_time_tuning_curves = [];
pooled_time_MI = [];

pooled_distance_pvalue = [];
pooled_distance_tuning_curves = [];
pooled_distance_MI = [];

%% Load zScored MI values
for mouse=1:length(mice)
    for session=1:length(session2use)
        cd(['/Volumes/GoogleDrive/My Drive/Williams lab/Manuscripts/ETTER2022_NatComs/input_data/calcium_imaging/' mice{mouse} '/LT/LongMonoTrack/' mice{mouse} '_' session2use{session}])
        if isfile('LMT_results.mat')
            load('LMT_results.mat');
            pooled_spatial_tuning_curves = [pooled_spatial_tuning_curves spatial_tuning_curve_data];
            pooled_spatial_pvalue = [pooled_spatial_pvalue spatial_MI_pvalue];
            pooled_spatial_MI = [pooled_spatial_MI spatial_MI];
            
            pooled_time_pvalue = [pooled_time_pvalue time_MI_pvalue];
            pooled_time_tuning_curves = [pooled_time_tuning_curves time_tuning_curve_data];
            pooled_time_MI = [pooled_time_MI time_MI];
            
            pooled_distance_pvalue = [pooled_distance_pvalue distance_MI_pvalue];
            pooled_distance_tuning_curves = [pooled_distance_tuning_curves distance_tuning_curve_data];
            pooled_distance_MI = [pooled_distance_MI distance_MI];
        end
    end
end

%% Smooth
for i=1:size(pooled_spatial_tuning_curves,2)
    smoothed_spatial_tuning_curves(:,i) = smoothPlaceField4Plots(pooled_spatial_tuning_curves(:,i),5);
    smoothed_time_tuning_curves(:,i) = smoothPlaceField4Plots(pooled_time_tuning_curves(:,i),10);
    smoothed_distance_tuning_curves(:,i) = smoothPlaceField4Plots(pooled_distance_tuning_curves(:,i),15);
end
    
%% Sort
[sorted_spatial_MI, sorted_space_MI_idx] = sort(pooled_spatial_MI, 'descend');
[sorted_time_MI, sorted_time_MI_idx] = sort(pooled_time_MI, 'descend');
[sorted_distance_MI, sorted_distance_MI_idx] = sort(pooled_distance_MI, 'descend');

% for neuron=1:length(pooled_spatial_MI)
%    plot(smoothed_spatial_tuning_curves(:,sorted_space_MI_idx(neuron)));
%    ax=gca;
%    ax.YLim = [0 0.3];
%    neuronStr = {'Neuron: ' neuron 'MI: ' sorted_spatial_MI(neuron) 'pvalue:' pooled_spatial_pvalue(sorted_space_MI_idx(neuron))};
%    title(neuronStr)
%    pause
% end

% Spatial neurons:
spatial_rank1=13; % MI=0.080, p-value=0.001
spatial_rank2=105; % MI=0.053, p-value=0.001
spatial_rank3=290; % MI=0.043, p-value=0.004
spatial_rank4=672; % MI=0.033, p-value=0.011
spatial_rank5=929; % MI=0.029, p-value=0.154
spatial_rank6=1259; % MI=0.026, p-value=0.511

% for neuron=1:length(pooled_time_MI)
%    plot(smoothed_time_tuning_curves(:,sorted_time_MI_idx(neuron)));
%    ax=gca;
%    ax.YLim = [0 0.2];
%    neuronStr = {'Neuron: ' neuron 'MI: ' sorted_time_MI(neuron) 'pvalue:' pooled_time_pvalue(sorted_time_MI_idx(neuron))};
%    title(neuronStr)
%    pause
% end

temporal_rank1=3; % MI=0.076, p-value<0.0001
temporal_rank2=61; % MI=0.053, p-value<0.0001
temporal_rank3=131; % MI=0.045, p-value=0.004
%temporal_rank3=132; % MI=0.045, p-value=0.023
temporal_rank4=188; % MI=0.042, p-value=0.012
%temporal_rank5bis=228; % MI=0.040, p-value=0.021
temporal_rank5=250; % MI=0.039, p-value=0.055
temporal_rank6=590; % MI=0.031, p-value=0.366

% for neuron=1:length(pooled_distance_MI)
%    plot(smoothed_distance_tuning_curves(:,sorted_distance_MI_idx(neuron)));
%    ax=gca;
%    ax.YLim = [0 1];
%    neuronStr = {'Neuron: ' neuron 'MI: ' sorted_distance_MI(neuron) 'pvalue:' pooled_distance_pvalue(sorted_distance_MI_idx(neuron))};
%    title(neuronStr)
%    pause
% end

distance_rank1=7; % MI=0.139, p-value=0.001
distance_rank2=44; % MI=0.109, p-value=0.008
distance_rank3=102; % MI=0.094, p-value=0.005
%distance_rank3bis=96; % MI=0.095, p-value=0.027
distance_rank4=112; % MI=0.092, p-value=0.001
distance_rank5=131; % MI=0.086, p-value=0.015
distance_rank6=423; % MI=0.064, p-value=0.164


%% Space
figure
subplot(7,1,1);
plot(sorted_spatial_MI); hold on;
line([spatial_rank1 spatial_rank1],[0 1])
line([spatial_rank2 spatial_rank2],[0 1])
line([spatial_rank3 spatial_rank3],[0 1])
line([spatial_rank4 spatial_rank4],[0 1])
line([spatial_rank5 spatial_rank5],[0 1])
line([spatial_rank6 spatial_rank6],[0 1])
ax=gca;
ax.XLim=[1 1500];
ax.YLim=[0 0.2];

subplot(7,1,2);
title('Spatial');
plot(smoothed_spatial_tuning_curves(:,sorted_space_MI_idx(spatial_rank1)));
ax=gca;
ax.YLim=[0 0.3];

subplot(7,1,3);
plot(smoothed_spatial_tuning_curves(:,sorted_space_MI_idx(spatial_rank2)));
ax=gca;
ax.YLim=[0 0.3];

subplot(7,1,4);
plot(smoothed_spatial_tuning_curves(:,sorted_space_MI_idx(spatial_rank3)));
ax=gca;
ax.YLim=[0 0.3];

subplot(7,1,5);
plot(smoothed_spatial_tuning_curves(:,sorted_space_MI_idx(spatial_rank4)));
ax=gca;
ax.YLim=[0 0.3];

subplot(7,1,6);
plot(smoothed_spatial_tuning_curves(:,sorted_space_MI_idx(spatial_rank5)));
ax=gca;
ax.YLim=[0 0.3];

subplot(7,1,7);
plot(smoothed_spatial_tuning_curves(:,sorted_space_MI_idx(spatial_rank6)));
ax=gca;
ax.YLim=[0 0.3];

%% Time
figure
subplot(7,1,1);
plot(sorted_time_MI); hold on;
line([temporal_rank1 temporal_rank1],[0 1])
line([temporal_rank2 temporal_rank2],[0 1])
line([temporal_rank3 temporal_rank3],[0 1])
line([temporal_rank4 temporal_rank4],[0 1])
line([temporal_rank5 temporal_rank5],[0 1])
line([temporal_rank6 temporal_rank6],[0 1])
ax=gca;
ax.XLim=[1 1000];
ax.YLim=[0 0.2];

subplot(7,1,2);
title('Spatial');
plot(smoothed_time_tuning_curves(:,sorted_time_MI_idx(temporal_rank1)));
ax=gca;
ax.XLim=[0 50];
ax.YLim=[0 0.2];

subplot(7,1,3);
plot(smoothed_time_tuning_curves(:,sorted_time_MI_idx(temporal_rank2)));
ax=gca;
ax.XLim=[0 50];
ax.YLim=[0 0.2];

subplot(7,1,4);
plot(smoothed_time_tuning_curves(:,sorted_time_MI_idx(temporal_rank3)));
ax=gca;
ax.XLim=[0 50];
ax.YLim=[0 0.2];

subplot(7,1,5);
plot(smoothed_time_tuning_curves(:,sorted_time_MI_idx(temporal_rank4)));
ax=gca;
ax.XLim=[0 50];
ax.YLim=[0 0.2];

subplot(7,1,6);
plot(smoothed_time_tuning_curves(:,sorted_time_MI_idx(temporal_rank5)));
ax=gca;
ax.XLim=[0 50];
ax.YLim=[0 0.2];

subplot(7,1,7);
plot(smoothed_time_tuning_curves(:,sorted_time_MI_idx(temporal_rank6)));
ax=gca;
ax.XLim=[0 50];
ax.YLim=[0 0.2];

%% Distance
figure
subplot(7,1,1);
plot(sorted_distance_MI); hold on;
line([distance_rank1 distance_rank1],[0 1])
line([distance_rank2 distance_rank2],[0 1])
line([distance_rank3 distance_rank3],[0 1])
line([distance_rank4 distance_rank4],[0 1])
line([distance_rank5 distance_rank5],[0 1])
line([distance_rank6 distance_rank6],[0 1])
ax=gca;
ax.XLim=[1 500];
ax.YLim=[0 0.2];

subplot(7,1,2);
title('Distance');
plot(smoothed_distance_tuning_curves(:,sorted_distance_MI_idx(distance_rank1)));
ax=gca;
ax.YLim=[0 0.8];

subplot(7,1,3);
plot(smoothed_distance_tuning_curves(:,sorted_distance_MI_idx(distance_rank2)));
ax=gca;
ax.YLim=[0 0.8];

subplot(7,1,4);
plot(smoothed_distance_tuning_curves(:,sorted_distance_MI_idx(distance_rank3)));
ax=gca;
ax.YLim=[0 0.8];

subplot(7,1,5);
plot(smoothed_distance_tuning_curves(:,sorted_distance_MI_idx(distance_rank4)));
ax=gca;
ax.YLim=[0 0.8];

subplot(7,1,6);
plot(smoothed_distance_tuning_curves(:,sorted_distance_MI_idx(distance_rank5)));
ax=gca;
ax.YLim=[0 0.8];

subplot(7,1,7);
plot(smoothed_distance_tuning_curves(:,sorted_distance_MI_idx(distance_rank6)));
ax=gca;
ax.YLim=[0 0.8];

end

