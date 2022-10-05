function [zMI, MI, p_value] = extractCrossMI(input_matrix, plotting)
%EXTRACTCROSSMI Summary of this function goes here
%   Takes a matrix of binarized activity, and spits out matrices of MI and
%   significance

numShuffles = 1000;
numFrames = size(input_matrix,1);
numNeurons = size(input_matrix,2);

zMI=zeros(size(input_matrix,2))*nan;
MI=zeros(size(input_matrix,2))*nan;
p_value=zeros(size(input_matrix,2))*nan;


for cell_i = 1:numNeurons
    disp(['Progress: ' num2str(cell_i/numNeurons*100) '%']);
    for cell_j = 1:numNeurons
        
        if cell_j <= cell_i % This is to avoid computing the symmetric part (and save time)
            % First compute marginal probabilities that will not change
            % after shuffling:
            trace_i = input_matrix(:,cell_i);
            trace_j = input_matrix(:,cell_j);
                
            [MI(cell_i,cell_j), p_value(cell_i,cell_j), zMI(cell_i,cell_j)] = compute_MI_significance(trace_i,trace_j);
            %MI(cell_j,cell_i) = MI(cell_i,cell_j); % For symmetry
            %p_value(cell_j,cell_i) = p_value(cell_i,cell_j); % For symmetry
            %zMI(cell_j,cell_i) = zMI(cell_i,cell_j); % For symmetry
            
            if plotting
            imagesc(zMI)
            set(gca, 'CLim',[2 10]);
            colorbar
            drawnow
            end
            
        end
    end
end
end

