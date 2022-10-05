function [velocity] = extract_DLC_velocity(input_structure)
%EXTRACT_DLC_VELOCITY Summary of this function goes here
%   Detailed explanation goes here

time_vec=0;
Fs = 30;
dt = 1/Fs;
position = input_structure.position;
for i = 2:size(position,1)
   time_vec(i) =  time_vec(i-1) + dt; 
end

velocity(1,:) = 0;
    for i=2:length(time_vec)
        velocity(i,:) = sqrt((position(i,1)-position(i-1,1)).^2 + (position(i,2)-position(i-1,2)).^2)/(time_vec(i)-time_vec(i-1));
    end
    
velocity = smooth(velocity,round(1/mode(diff(time_vec))));

end

