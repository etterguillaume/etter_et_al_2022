function [MI, p_value, zMI] = compute_MI_significance(vec1, vec2)
    % vec1 and vec2 should be binary vectors of equal length

vecLength = length(vec1);           
            
% First, extract the actual MI from your binarized activity and inclusion vector (vec2)              
MI = extractBinaryMI(vec1,vec2);

% Then do your shuffling here (1000x at least to get a precise p-value)
numShuffles = 1000;
shuffled_MI = zeros(numShuffles,1); % Initialize matrix value
for shuffle_i = 1:numShuffles
    random_ts = ceil(rand*vecLength);
    shuffledVec1 = zeros(vecLength,1);
    shuffledVec1(1:random_ts) = vec1(end-random_ts+1:end);
    shuffledVec1(random_ts+1:end) = vec1(1:end-random_ts);
    
    shuffled_MI(shuffle_i) = extractBinaryMI(shuffledVec1,vec2);
end

p_value = sum(shuffled_MI > MI)/numShuffles;

zMI = (MI-mean(shuffled_MI))/std(shuffled_MI);