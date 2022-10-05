function HDist = extract_CrossHDist(input_matrix, plotting)
%EXTRACTCROSSMI Summary of this function goes here
%   Takes a matrix of binarized activity, and spits out matrices of MI and
%   significance

numNeurons = size(input_matrix,2);

HDist=zeros(size(input_matrix,2))*nan;

for cell_i = 1:numNeurons
    disp(['Progress: ' num2str(cell_i/numNeurons*100) '%']);
    for cell_j = 1:numNeurons
        
        if cell_j < cell_i % This is to avoid computing the symmetric part (and save time)
            % First compute marginal probabilities that will not change
            % after shuffling:
            trace_i = input_matrix(:,cell_i);
            trace_j = input_matrix(:,cell_j);
                
            HDist(cell_i,cell_j) = extract_HammingDistance(trace_i,trace_j);
            
            if plotting
            imagesc(HDist)
            daspect([1 1 1])
            colorbar
            colormap Viridis
            drawnow
            end
            
        end
    end
end
end

