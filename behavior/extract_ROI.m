function ROIs = extract_ROI(behav)
%EXTRACT_ROI Summary of this function goes here
%   Detailed explanation goes here

ROIs={};

figure
imshow(behav.background)

user_resp = 'y'

while user_resp == 'y'
current_ROI_name = input(['Name this ROI: '], 's');

disp('Define this ROI then press Return');
current_ROI = drawrectangle('label', current_ROI_name, 'Color',[rand rand rand]);

current_ROI_dims = current_ROI.Position;

ROIs(end+1).name = current_ROI_name;
ROIs(end).dims = current_ROI_dims;

test = 0;

user_resp = input(['Do you want to define another ROIs? (y/n): '], 's');

end

close all

end

