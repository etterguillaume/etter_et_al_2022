function [MI] = extractBinaryMI(binary1,binary2)
%EXTRACTBINARYMI Summary of this function goes here
%   Extract mutual information between two binary vectors of same length
%   Guillaume Etter

binary1 = logical(binary1);
binary2 = logical(binary2);

if length(binary1) ~= length(binary2)
    error('Input vectors must be the same length')
else 
lengthBinary = length(binary1);
end

MI = 0; % Initialize matrix value
for state_i = 0:1
    for state_j = 0:1
        pState_i = sum(binary1 == state_i)/lengthBinary;
        pState_j = sum(binary2 == state_j)/lengthBinary;
        jointProb = sum(binary1 == state_i & binary2 == state_j)/lengthBinary;
        if jointProb ~= 0
            MI = MI + jointProb*log2(jointProb./(pState_i*pState_j));
        end
    end
end


end

