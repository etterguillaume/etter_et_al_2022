function plot_CellTypes
%PLOT_CELLTYPES Summary of this function goes here
%   Detailed explanation goes here

mice = ['PV986', 'PV988', 'PV989', 'PV990', 'PV991'];
session2use = 'LMT6_20190319';

%% Load zScored MI values
load(['/Volumes/GoogleDrive/My Drive/Williams lab/Manuscripts/ETTER2022_NatComs/input_data/calcium_imaging/PV986/LT/LongMonoTrack/PV986_' session2use '/LMT_results.mat']);
pooled_distance_pvalue = distance_MI_pvalue;
pooled_spatial_pvalue = spatial_MI_pvalue;
pooled_time_pvalue = time_MI_pvalue;
pooled_distance_tuning_curves = distance_tuning_curve_data;
pooled_spatial_tuning_curves = spatial_tuning_curve_data;
pooled_time_tuning_curves = time_tuning_curve_data;

load(['/Volumes/GoogleDrive/My Drive/Williams lab/Manuscripts/ETTER2022_NatComs/input_data/calcium_imaging/PV988/LT/LongMonoTrack/PV988_' session2use '/LMT_results.mat']);
pooled_distance_pvalue = [pooled_distance_pvalue distance_MI_pvalue];
pooled_spatial_pvalue = [pooled_spatial_pvalue spatial_MI_pvalue];
pooled_time_pvalue = [pooled_time_pvalue time_MI_pvalue];
pooled_distance_tuning_curves = [pooled_distance_tuning_curves distance_tuning_curve_data];
pooled_spatial_tuning_curves = [pooled_spatial_tuning_curves spatial_tuning_curve_data];
pooled_time_tuning_curves = [pooled_time_tuning_curves time_tuning_curve_data];

load(['/Volumes/GoogleDrive/My Drive/Williams lab/Manuscripts/ETTER2022_NatComs/input_data/calcium_imaging/PV989/LT/LongMonoTrack/PV989_' session2use '/LMT_results.mat']);
pooled_distance_pvalue = [pooled_distance_pvalue distance_MI_pvalue];
pooled_spatial_pvalue = [pooled_spatial_pvalue spatial_MI_pvalue];
pooled_time_pvalue = [pooled_time_pvalue time_MI_pvalue];
pooled_distance_tuning_curves = [pooled_distance_tuning_curves distance_tuning_curve_data];
pooled_spatial_tuning_curves = [pooled_spatial_tuning_curves spatial_tuning_curve_data];
pooled_time_tuning_curves = [pooled_time_tuning_curves time_tuning_curve_data];

%load(['/Volumes/GoogleDrive/My Drive/Williams lab/Williams lab/Manuscripts/ETTER2022_NatComs/input_data/calcium_imaging/PV990/LT/LongMonoTrack/PV990_' session2use '/LMT_results.mat']);
%pooled_distance_pvalue = [pooled_distance_pvalue distance_MI_pvalue];
%pooled_spatial_pvalue = [pooled_spatial_pvalue spatial_MI_pvalue];
%pooled_time_pvalue = [pooled_time_pvalue time_MI_pvalue];
%pooled_distance_tuning_curves = [pooled_distance_tuning_curves
%distance_tuning_curve_data];
%pooled_spatial_tuning_curves = [pooled_spatial_tuning_curves spatial_tuning_curve_data];
%pooled_time_tuning_curves = [pooled_time_tuning_curves time_tuning_curve_data];

load(['/Volumes/GoogleDrive/My Drive/Williams lab/Manuscripts/ETTER2022_NatComs/input_data/calcium_imaging/PV991/LT/LongMonoTrack/PV991_' session2use '/LMT_results.mat']);
pooled_distance_pvalue = [pooled_distance_pvalue distance_MI_pvalue];
pooled_spatial_pvalue = [pooled_spatial_pvalue spatial_MI_pvalue];
pooled_time_pvalue = [pooled_time_pvalue time_MI_pvalue];
pooled_distance_tuning_curves = [pooled_distance_tuning_curves distance_tuning_curve_data];
pooled_spatial_tuning_curves = [pooled_spatial_tuning_curves spatial_tuning_curve_data];
pooled_time_tuning_curves = [pooled_time_tuning_curves time_tuning_curve_data];

figure
%% Obtain sorting indices
%% Place cells
nonPlaceCells = pooled_spatial_pvalue >= .05 | pooled_time_pvalue<0.05 | pooled_distance_pvalue<0.05;
sigSpace_tuning_curve_data = pooled_spatial_tuning_curves;
sigSpace_time_tCurve = pooled_time_tuning_curves;
sigSpace_distance_tCurve = pooled_distance_tuning_curves;

sigSpace_tuning_curve_data(:,nonPlaceCells) = []; % Remove non place cells
sigSpace_time_tCurve(:,nonPlaceCells) = []; % Remove non place cells
sigSpace_distance_tCurve(:,nonPlaceCells) = []; % Remove non place cells

%% Smooth
for i=1:size(sigSpace_distance_tCurve,2)
    sigSpace_distance_tCurve(:,i) = smoothPlaceFieldParam(sigSpace_distance_tCurve(:,i),15);
    sigSpace_time_tCurve(:,i) = smoothPlaceFieldParam(sigSpace_time_tCurve(:,i),10);
    sigSpace_tuning_curve_data(:,i) = smoothPlaceFieldParam(sigSpace_tuning_curve_data(:,i),5);
end
    
%% Sort
[~, peak_space_index] = max(sigSpace_tuning_curve_data,[],1);
[~, sorted_space_index] = sort(peak_space_index);

% Plot place cell data, row 1
subplot(3,3,1)
imagesc(sigSpace_tuning_curve_data(:,sorted_space_index)');
ax1=gca;
ax1.CLim = [0 0.15];
daspect([1 1 1])
colormap Viridis
colorbar
title 'Space'

subplot(3,3,2)
imagesc(sigSpace_time_tCurve(:,sorted_space_index)');
ax2=gca;
ax2.CLim = [0 0.15];
daspect([1 1 1])
colormap Viridis
colorbar
title 'Time'

subplot(3,3,3)
imagesc(sigSpace_distance_tCurve(:,sorted_space_index)');
ax3=gca;
ax3.CLim = [0 0.15];
daspect([1 1 1])
colormap Viridis
colorbar
title 'Distance'

%% Time cells
nonTimeCells = pooled_time_pvalue >= 0.05 | pooled_spatial_pvalue<0.05 | pooled_distance_pvalue<0.05;
sigTime_tuning_curve_data = pooled_time_tuning_curves;
sigTime_place_tCurve = pooled_spatial_tuning_curves;
sigTime_distance_tCurve = pooled_distance_tuning_curves;

sigTime_tuning_curve_data(:,nonTimeCells) = []; % Remove non place cells
sigTime_place_tCurve(:,nonTimeCells) = []; % Remove non place cells
sigTime_distance_tCurve(:,nonTimeCells) = []; % Remove non place cells

%% Smooth
for i=1:size(sigTime_distance_tCurve,2)
    sigTime_distance_tCurve(:,i) = smoothPlaceField(sigTime_distance_tCurve(:,i));
    sigTime_place_tCurve(:,i) = smoothPlaceField(sigTime_place_tCurve(:,i));
    sigTime_tuning_curve_data(:,i) = smoothPlaceField(sigTime_tuning_curve_data(:,i));
end

%% Sort
[~, peak_time_index] = max(sigTime_tuning_curve_data,[],1);
[~, sorted_time_index] = sort(peak_time_index);

%% Plot place cell data, row 1
subplot(3,3,4)
imagesc(sigTime_place_tCurve(:,sorted_time_index)');
ax4=gca;
ax4.CLim = [0 0.2];
daspect([1 1 1])
colormap Viridis
colorbar
title 'Space'

subplot(3,3,5)
imagesc(sigTime_tuning_curve_data(:,sorted_time_index)');
ax5=gca;
ax5.CLim = [0 0.2];
daspect([1 1 1])
colormap Viridis
colorbar
title 'Time'

subplot(3,3,6)
imagesc(sigTime_distance_tCurve(:,sorted_time_index)');
ax6=gca;
ax6.CLim = [0 0.2];
daspect([1 1 1])
colormap Viridis
colorbar
title 'Distance'


%% Distance cells
nonDistanceCells = pooled_distance_pvalue >= 0.05 | pooled_spatial_pvalue<0.05 | pooled_time_pvalue<0.05;
sigDistance_tuning_curve_data = pooled_distance_tuning_curves;
sigDistance_place_tCurve = pooled_spatial_tuning_curves;
sigDistance_time_tCurve = pooled_time_tuning_curves;

sigDistance_tuning_curve_data(:,nonDistanceCells) = []; % Remove non place cells
sigDistance_place_tCurve(:,nonDistanceCells) = []; % Remove non place cells
sigDistance_time_tCurve(:,nonDistanceCells) = []; % Remove non place cells

%% Smooth
for i=1:size(sigDistance_tuning_curve_data,2)
    sigDistance_place_tCurve(:,i) = smoothPlaceField(sigDistance_place_tCurve(:,i));
    sigDistance_time_tCurve(:,i) = smoothPlaceField(sigDistance_time_tCurve(:,i));
    sigDistance_tuning_curve_data(:,i) = smoothPlaceField(sigDistance_tuning_curve_data(:,i));
end

%% Sort
[~, peak_distance_index] = max(sigDistance_tuning_curve_data,[],1);
[~, sorted_distance_index] = sort(peak_distance_index);

% Plot place cell data, row 1
subplot(3,3,7)
imagesc(sigDistance_place_tCurve(:,sorted_distance_index)');
ax7=gca;
ax7.CLim = [0 0.2];
daspect([1 1 1])
colormap Viridis
colorbar
title 'Space'

subplot(3,3,8)
imagesc(sigDistance_time_tCurve(:,sorted_distance_index)');
ax8=gca;
ax8.CLim = [0 0.2];
daspect([1 1 1])
colormap Viridis
colorbar
title 'Time'

subplot(3,3,9)
imagesc(sigDistance_tuning_curve_data(:,sorted_distance_index)');
ax9=gca;
ax9.CLim = [0 0.2];
daspect([1 1 1])
colormap Viridis
colorbar
title 'Distance'

end

