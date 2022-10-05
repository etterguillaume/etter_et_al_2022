workingDir = pwd;
load([workingDir '/LMT_results.mat'])

placeCellsProper = spatial_MI_pvalue < 0.05 & time_MI_pvalue > 0.05 & distance_MI_pvalue>0.05;
timeCellsProper = time_MI_pvalue < 0.05 & spatial_MI_pvalue>0.05 & distance_MI_pvalue>0.05;
distanceCellsProper = distance_MI_pvalue < 0.05 & spatial_MI_pvalue > 0.05 & time_MI_pvalue > 0.05;

conjSpatial = spatial_MI_pvalue < 0.05;
conjTemporal = time_MI_pvalue < 0.05;
conjDistance = distance_MI_pvalue < 0.05;

nonSpatial = spatial_MI_pvalue > 0.05;
nonTemporal = time_MI_pvalue > 0.05;
nonDistance = distance_MI_pvalue > 0.05;

disp(['Place cells: ' num2str(sum(placeCellsProper))])
disp(['Time cells: ' num2str(sum(timeCellsProper))])
disp(['Distance cells: ' num2str(sum(distanceCellsProper))])

disp(['Conjunctive spatial cells: ' num2str(sum(conjSpatial))])
disp(['Conjunctive temporal cells: ' num2str(sum(conjTemporal))])
disp(['Conjunctive distance cells: ' num2str(sum(conjDistance))])

disp(['Non spatial cells: ' num2str(sum(nonSpatial))])
disp(['Non temporal cells: ' num2str(sum(nonTemporal))])
disp(['Non distance cells: ' num2str(sum(nonDistance))])

