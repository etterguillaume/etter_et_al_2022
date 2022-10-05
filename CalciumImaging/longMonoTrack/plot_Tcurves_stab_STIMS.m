function plot_CellTypes
%PLOT_CELLTYPES Summary of this function goes here
%   Detailed explanation goes here

%% Load zScored MI values
% PV986
load(['/Volumes/GoogleDrive/My Drive/Williams lab/Manuscripts/ETTER_NatNeurosc/input_data/calcium_imaging/PV986/LT/LongMonoTrack/LMT_spatial_stability_12_35_46_PCs_only.mat']);
pooled_PF1 = placeFields1;
pooled_PF2 = placeFields2;
pooled_PF3 = placeFields3;
pooled_PF4 = placeFields4;
pooled_PF5 = placeFields5;
pooled_PF6 = placeFields6;
pooled_PF_stab_12 = stability_12;
pooled_PF_stab_35 = stability_35;
pooled_PF_stab_46 = stability_46;
load(['/Volumes/GoogleDrive/My Drive/Williams lab/Manuscripts/ETTER_NatNeurosc/input_data/calcium_imaging/PV986/LT/LongMonoTrack/LMT_temporal_stability_12_35_46_TCs_only.mat']);
pooled_TF1 = timeFields1;
pooled_TF2 = timeFields2;
pooled_TF3 = timeFields3;
pooled_TF4 = timeFields4;
pooled_TF5 = timeFields5;
pooled_TF6 = timeFields6;
pooled_TF_stab_12 = stability_12;
pooled_TF_stab_35 = stability_35;
pooled_TF_stab_46 = stability_46;

% PV989
%load(['/Users/guillaumeetter/Google Drive/Williams lab/Manuscripts/ETTER_NatNeurosc/input_data/calcium_imaging/PV989/LT/LongMonoTrack/LMT_spatial_stability_12_35_46_PCs_only.mat']);
% pooled_PF1 = [pooled_PF1 placeFields1];
% pooled_PF2 = [pooled_PF2 placeFields2];
% pooled_PF3 = [pooled_PF3 placeFields3];
% pooled_PF4 = [pooled_PF4 placeFields4];
% pooled_PF5 = [pooled_PF5 placeFields5];
% pooled_PF6 = [pooled_PF6 placeFields6];
% pooled_PF_stab_12 = [pooled_PF_stab_12 stability_12];
% pooled_PF_stab_35 = [pooled_PF_stab_35 stability_35];
% pooled_PF_stab_46 = [pooled_PF_stab_46 stability_46];
% load(['/Users/guillaumeetter/Google Drive/Williams lab/Manuscripts/ETTER_NatNeurosc/input_data/calcium_imaging/PV989/LT/LongMonoTrack/LMT_temporal_stability_12_35_46_TCs_only.mat']);
% pooled_TF1 = [pooled_TF1 timeFields1];
% pooled_TF2 = [pooled_TF2 timeFields2];
% pooled_TF3 = [pooled_TF3 timeFields3];
% pooled_TF4 = [pooled_TF4 timeFields4];
% pooled_TF5 = [pooled_TF5 timeFields5];
% pooled_TF6 = [pooled_TF6 timeFields6];
% pooled_TF_stab_12 = [pooled_TF_stab_12 stability_12];
% pooled_TF_stab_35 = [pooled_TF_stab_35 stability_35];
% pooled_TF_stab_46 = [pooled_TF_stab_46 stability_46];

% PV990
load(['/Volumes/GoogleDrive/My Drive/Williams lab/Manuscripts/ETTER_NatNeurosc/input_data/calcium_imaging/PV990/LT/LongMonoTrack/LMT_spatial_stability_12_35_46_PCs_only.mat']);
pooled_PF1 = [pooled_PF1 placeFields1];
pooled_PF2 = [pooled_PF2 placeFields2];
pooled_PF3 = [pooled_PF3 placeFields3];
pooled_PF4 = [pooled_PF4 placeFields4];
pooled_PF5 = [pooled_PF5 placeFields5];
pooled_PF6 = [pooled_PF6 placeFields6];
pooled_PF_stab_12 = [pooled_PF_stab_12 stability_12];
pooled_PF_stab_35 = [pooled_PF_stab_35 stability_35];
pooled_PF_stab_46 = [pooled_PF_stab_46 stability_46];
load(['/Volumes/GoogleDrive/My Drive/Williams lab/Manuscripts/ETTER_NatNeurosc/input_data/calcium_imaging/PV990/LT/LongMonoTrack/LMT_temporal_stability_12_35_46_TCs_only.mat']);
pooled_TF1 = [pooled_TF1 timeFields1];
pooled_TF2 = [pooled_TF2 timeFields2];
pooled_TF3 = [pooled_TF3 timeFields3];
pooled_TF4 = [pooled_TF4 timeFields4];
pooled_TF5 = [pooled_TF5 timeFields5];
pooled_TF6 = [pooled_TF6 timeFields6];
pooled_TF_stab_12 = [pooled_TF_stab_12 stability_12];
pooled_TF_stab_35 = [pooled_TF_stab_35 stability_35];
pooled_TF_stab_46 = [pooled_TF_stab_46 stability_46];

% PV991
load(['/Volumes/GoogleDrive/My Drive/Williams lab/Manuscripts/ETTER_NatNeurosc/input_data/calcium_imaging/PV991/LT/LongMonoTrack/LMT_spatial_stability_12_35_46_PCs_only.mat']);
pooled_PF1 = [pooled_PF1 placeFields1];
pooled_PF2 = [pooled_PF2 placeFields2];
pooled_PF3 = [pooled_PF3 placeFields3];
pooled_PF4 = [pooled_PF4 placeFields4];
pooled_PF5 = [pooled_PF5 placeFields5];
pooled_PF6 = [pooled_PF6 placeFields6];
pooled_PF_stab_12 = [pooled_PF_stab_12 stability_12];
pooled_PF_stab_35 = [pooled_PF_stab_35 stability_35];
pooled_PF_stab_46 = [pooled_PF_stab_46 stability_46];
load(['/Volumes/GoogleDrive/My Drive/Williams lab/Manuscripts/ETTER_NatNeurosc/input_data/calcium_imaging/PV991/LT/LongMonoTrack/LMT_temporal_stability_12_35_46_TCs_only.mat']);
pooled_TF1 = [pooled_TF1 timeFields1];
pooled_TF2 = [pooled_TF2 timeFields2];
pooled_TF3 = [pooled_TF3 timeFields3];
pooled_TF4 = [pooled_TF4 timeFields4];
pooled_TF5 = [pooled_TF5 timeFields5];
pooled_TF6 = [pooled_TF6 timeFields6];
pooled_TF_stab_12 = [pooled_TF_stab_12 stability_12];
pooled_TF_stab_35 = [pooled_TF_stab_35 stability_35];
pooled_TF_stab_46 = [pooled_TF_stab_46 stability_46];

%% Compute probability density function
% optional
% Note that here tuning curves should already be smoothed and pc/tc
% identified and selected
    
%% Sort
[~, peak_PF1_index] = max(pooled_PF1,[],1);
[~, sorted_PF1_index] = sort(peak_PF1_index);
[~, peak_PF3_index] = max(pooled_PF3,[],1);
[~, sorted_PF3_index] = sort(peak_PF3_index);
[~, peak_PF4_index] = max(pooled_PF4,[],1);
[~, sorted_PF4_index] = sort(peak_PF4_index);

[~, peak_TF1_index] = max(pooled_TF1,[],1);
[~, sorted_TF1_index] = sort(peak_TF1_index);
[~, peak_TF3_index] = max(pooled_TF3,[],1);
[~, sorted_TF3_index] = sort(peak_TF3_index);
[~, peak_TF4_index] = max(pooled_TF4,[],1);
[~, sorted_TF4_index] = sort(peak_TF4_index);

figure
subplot(3,3,1)
imagesc(pooled_PF1(:,sorted_PF1_index)');
ax1=gca;
%ax1.CLim = [0 0.2];
daspect([1 1 1])
colormap Viridis
colorbar
title 'PF - BL'

subplot(3,3,2)
imagesc(pooled_PF2(:,sorted_PF1_index)');
ax2=gca;
%ax2.CLim = [0 0.2];
daspect([1 1 1])
colormap Viridis
colorbar
title 'PF - no stim'

subplot(3,3,3)
plot(pooled_PF_stab_12(:,sorted_PF1_index));
ax3=gca;
view([90 -90])
title 'Stability'
%linkaxes([ax1 ax2 ax3], 'y')

subplot(3,3,4)
imagesc(pooled_PF3(:,sorted_PF3_index)');
ax4=gca;
%ax1.CLim = [0 0.2];
daspect([1 1 1])
colormap Viridis
colorbar
title 'PF - BL'

subplot(3,3,5)
imagesc(pooled_PF5(:,sorted_PF3_index)');
ax5=gca;
%ax2.CLim = [0 0.2];
daspect([1 1 1])
colormap Viridis
colorbar
title 'PF - 8Hz'

subplot(3,3,6)
plot(pooled_PF_stab_35(:,sorted_PF3_index));
ax6=gca;
view([90 -90])
title 'Stability'
%linkaxes([ax1 ax2 ax3], 'y')

subplot(3,3,7)
imagesc(pooled_PF4(:,sorted_PF4_index)');
ax7=gca;
%ax1.CLim = [0 0.2];
daspect([1 1 1])
colormap Viridis
colorbar
title 'PF - BL'

subplot(3,3,8)
imagesc(pooled_PF6(:,sorted_PF4_index)');
ax8=gca;
%ax2.CLim = [0 0.2];
daspect([1 1 1])
colormap Viridis
colorbar
title 'PF - Scrambled'

subplot(3,3,9)
plot(pooled_PF_stab_46(:,sorted_PF4_index)');
ax9=gca;
view([90 -90])
title 'Stability'

linkaxes([ax3 ax6 ax9], 'y')

%% Time cells
figure
subplot(3,3,1)
imagesc(pooled_TF1(:,sorted_TF1_index)');
ax1=gca;
%ax1.CLim = [0 0.2];
daspect([1 1 1])
colormap Viridis
colorbar
title 'TF - BL'

subplot(3,3,2)
imagesc(pooled_TF2(:,sorted_TF1_index)');
ax2=gca;
%ax2.CLim = [0 0.2];
daspect([1 1 1])
colormap Viridis
colorbar
title 'TF - no stim'

subplot(3,3,3)
plot(pooled_TF_stab_12(:,sorted_TF1_index));
ax3=gca;
view([90 -90])
title 'Stability'
%linkaxes([ax1 ax2 ax3], 'y')

subplot(3,3,4)
imagesc(pooled_TF3(:,sorted_TF3_index)');
ax4=gca;
%ax1.CLim = [0 0.2];
daspect([1 1 1])
colormap Viridis
colorbar
title 'TF - BL'

subplot(3,3,5)
imagesc(pooled_TF5(:,sorted_TF3_index)');
ax5=gca;
%ax2.CLim = [0 0.2];
daspect([1 1 1])
colormap Viridis
colorbar
title 'TF - 8Hz'

subplot(3,3,6)
plot(pooled_TF_stab_35(:,sorted_TF3_index));
ax6=gca;
view([90 -90])
title 'Stability'
%linkaxes([ax1 ax2 ax3], 'y')

subplot(3,3,7)
imagesc(pooled_TF4(:,sorted_TF4_index)');
ax7=gca;
%ax1.CLim = [0 0.2];
daspect([1 1 1])
colormap Viridis
colorbar
title 'TF - BL'

subplot(3,3,8)
imagesc(pooled_TF6(:,sorted_TF4_index)');
ax8=gca;
%ax2.CLim = [0 0.2];
daspect([1 1 1])
colormap Viridis
colorbar
title 'TF - Scrambled'

subplot(3,3,9)
plot(pooled_TF_stab_46(:,sorted_TF4_index)');
ax9=gca;
view([90 -90])
title 'Stability'

linkaxes([ax3 ax6 ax9], 'y')

end

