function plot_CellTypes_3dscatter
%PLOT_CELLTYPES_3DSCATTER Summary of this function goes here
%   Detailed explanation goes here

session2use = 'LMT6_20190319'

%% Load zScored MI values
load(['/Users/guillaumeetter/Google Drive/Williams lab/Manuscripts/ETTER_NatNeurosc/input_data/calcium_imaging/PV986/LT/LongMonoTrack/PV986_' session2use '/LMT_results.mat']);
pooled_distance_MI = distance_MI;
pooled_spatial_MI = spatial_MI;
pooled_time_MI = time_MI;
pooled_distance_pvalue = distance_MI_pvalue;
pooled_spatial_pvalue = spatial_MI_pvalue;
pooled_time_pvalue = time_MI_pvalue;
pooled_distance_zMI = distance_zMI;
pooled_spatial_zMI = spatial_zMI;
pooled_time_zMI = time_zMI;

load(['/Users/guillaumeetter/Google Drive/Williams lab/Manuscripts/ETTER_NatNeurosc/input_data/calcium_imaging/PV988/LT/LongMonoTrack/PV988_' session2use '/LMT_results.mat']);
pooled_distance_MI = [pooled_distance_MI distance_MI];
pooled_spatial_MI = [pooled_spatial_MI spatial_MI];
pooled_time_MI = [pooled_time_MI time_MI];
pooled_distance_pvalue = [pooled_distance_pvalue distance_MI_pvalue];
pooled_spatial_pvalue = [pooled_spatial_pvalue spatial_MI_pvalue];
pooled_time_pvalue = [pooled_time_pvalue time_MI_pvalue];
pooled_distance_zMI = [pooled_distance_zMI distance_zMI];
pooled_spatial_zMI = [pooled_spatial_zMI spatial_zMI];
pooled_time_zMI = [pooled_time_zMI time_zMI];

load(['/Users/guillaumeetter/Google Drive/Williams lab/Manuscripts/ETTER_NatNeurosc/input_data/calcium_imaging/PV989/LT/LongMonoTrack/PV989_' session2use '/LMT_results.mat']);
pooled_distance_MI = [pooled_distance_MI distance_MI];
pooled_spatial_MI = [pooled_spatial_MI spatial_MI];
pooled_time_MI = [pooled_time_MI time_MI];
pooled_distance_pvalue = [pooled_distance_pvalue distance_MI_pvalue];
pooled_spatial_pvalue = [pooled_spatial_pvalue spatial_MI_pvalue];
pooled_time_pvalue = [pooled_time_pvalue time_MI_pvalue];
pooled_distance_zMI = [pooled_distance_zMI distance_zMI];
pooled_spatial_zMI = [pooled_spatial_zMI spatial_zMI];
pooled_time_zMI = [pooled_time_zMI time_zMI];

load(['/Users/guillaumeetter/Google Drive/Williams lab/Manuscripts/ETTER_NatNeurosc/input_data/calcium_imaging/PV990/LT/LongMonoTrack/PV990_' session2use '/LMT_results.mat']);
pooled_distance_MI = [pooled_distance_MI distance_MI];
pooled_spatial_MI = [pooled_spatial_MI spatial_MI];
pooled_time_MI = [pooled_time_MI time_MI];
pooled_distance_pvalue = [pooled_distance_pvalue distance_MI_pvalue];
pooled_spatial_pvalue = [pooled_spatial_pvalue spatial_MI_pvalue];
pooled_time_pvalue = [pooled_time_pvalue time_MI_pvalue];
pooled_distance_zMI = [pooled_distance_zMI distance_zMI];
pooled_spatial_zMI = [pooled_spatial_zMI spatial_zMI];
pooled_time_zMI = [pooled_time_zMI time_zMI];

load(['/Users/guillaumeetter/Google Drive/Williams lab/Manuscripts/ETTER_NatNeurosc/input_data/calcium_imaging/PV991/LT/LongMonoTrack/PV991_' session2use '/LMT_results.mat']);
pooled_distance_MI = [pooled_distance_MI distance_MI];
pooled_spatial_MI = [pooled_spatial_MI spatial_MI];
pooled_time_MI = [pooled_time_MI time_MI];
pooled_distance_pvalue = [pooled_distance_pvalue distance_MI_pvalue];
pooled_spatial_pvalue = [pooled_spatial_pvalue spatial_MI_pvalue];
pooled_time_pvalue = [pooled_time_pvalue time_MI_pvalue];
pooled_distance_zMI = [pooled_distance_zMI distance_zMI];
pooled_spatial_zMI = [pooled_spatial_zMI spatial_zMI];
pooled_time_zMI = [pooled_time_zMI time_zMI];

%% Identify place/time cells proper
placeCellsProper = pooled_spatial_pvalue < 0.05 & pooled_time_pvalue > 0.05 & pooled_distance_pvalue>0.05;
timeCellsProper = pooled_time_pvalue < 0.05 & pooled_spatial_pvalue>0.05 & pooled_distance_pvalue>0.05;
distanceCellsProper = pooled_distance_pvalue < 0.05 & pooled_spatial_pvalue > 0.05 & pooled_time_pvalue > 0.05;

%sizeVec = placeCellsProper+timeCellsProper+distanceCellsProper;
%sizeVec2 = placeCellsProper+timeCellsProper;

figure
scatter(pooled_spatial_MI, pooled_spatial_zMI,10+placeCellsProper*20, [placeCellsProper' 0.5*placeCellsProper' 0.5*placeCellsProper']);
title 'Place cells'

figure
scatter(pooled_time_MI, pooled_time_zMI,10+timeCellsProper*20, [timeCellsProper' 0.5*timeCellsProper' 0.5*timeCellsProper']);
title 'Time cells'

figure
scatter(pooled_distance_MI, pooled_distance_zMI,10+distanceCellsProper*20, [distanceCellsProper' 0.5*distanceCellsProper' 0.5*distanceCellsProper']);
title 'Distance cells'
end

