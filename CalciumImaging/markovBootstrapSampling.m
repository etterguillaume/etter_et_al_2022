function [bootstrapTs] = markovBootstrapSampling(OFF2ON, ON2ON, numActiveLocs, inputTsLength)
%MARKOVBOOTSTRAPSAMPLING Creates a vector of random bootstrap samples using
%previously calculated transition probabilities
% Inputs:
% OFF2ON: probability of transitioning from 0 to 1
% ON2ON: probability of transitioning from 1 to 1
% numActiveLocs: total number of frames to include
% inputTsLength: length of the final bootstrap inclusion vector
% Outputs:
% bootstrapTs: bootstrap sampling vector (logical)

bootstrapTs = zeros(1,inputTsLength);

currentFrame = ceil(rand*inputTsLength); % Initialize a random frame
previousState = 0;

iter = 0;

while iter < inputTsLength || sum(bootstrapTs) < numActiveLocs
    iter = iter + 1;
    
    if previousState == 0
    currentState = rand < OFF2ON;
    elseif previousState == 1
    currentState = rand < ON2ON;
    end
    
    
    bootstrapTs(currentFrame) = currentState;
    %if currentState == 0
    %    currentFrame = ceil(rand*inputTsLength);
    %end
    
    previousState = currentState;
    
    currentFrame = currentFrame + 1;
    if currentFrame > inputTsLength
       currentFrame = 1; 
    end
    
    plot(bootstrapTs);
    drawnow
    
end

end

