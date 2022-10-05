function [MI, p_values, zMI, state_dependency] = MI_state_dependency_analysis(activity_matrix, states_matrix, states_names)
%MI_STATE_ANALYSIS Summary of this function goes here
%   INPUTS:
%   activity_matrix: t x n matrix of binarized activity over time t for each
%   neuron
%   states_matrix: a t x s matrix of states s over time t
%   OUTPUTS:
%   MI: raw MI values
%   p_values: significance against circularly permuted data
%   z_values: z-scored MI using circularly permuted data
%   state_dependency: a logical matrix indicating if states are dependent (1) or
%   not (0)

numFrames = size(activity_matrix,1);
numStates = size(states_matrix,2);
state_dependency = zeros(numStates);

for state_i = 1:numStates
    for state_j = 1:numStates
        if state_j > state_i % To avoid redundancy
            if sum(states_matrix(:,state_i) == 1 & states_matrix(:,state_j) == 1)/numFrames == 0
            state_dependency(state_i,state_j) = 0;
            else
            state_dependency(state_i,state_j) = 1;
            end
            state_dependency(state_j,state_i) = state_dependency(state_i,state_j); % For symmetry
        end
    end
end

state_list = [];
for state_interaction_level = 1:numStates
    for state_i = 1:state_interaction_level*numStates
    state_list(end+1) = state_names(state_i);
    end
    
    

end


end

