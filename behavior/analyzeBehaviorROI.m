function [results] = analyzeBehaviorROI(behav, ROIs)
%ANALYSEBEHAVIOR Summary of this function goes here
%   Detailed explanation goes here

numROIs = length(ROIs);
position = behav.nosePosition;
position = fillmissing(position,'linear', 'EndValues', 'nearest');

time_vec = behav.time/1000;

delta_t = mode(diff(time_vec));

for ROI_i = 1:numROIs
   current_ROI_name = ROIs(ROI_i).name;
   current_ROI_dims = ROIs(ROI_i).dims*behav.cmPerPixels;
   
   results.(current_ROI_name).inROI_index = find(position(:,1) > current_ROI_dims(1) & position(:,2) > current_ROI_dims(2) & position(:,1) < current_ROI_dims(1)+current_ROI_dims(3) & position(:,2) < current_ROI_dims(2)+current_ROI_dims(4));
   results.(current_ROI_name).inROI_time = delta_t*length(results.(current_ROI_name).inROI_index);
   
    
end

end

