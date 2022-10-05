function [MI, p_values, zMI] = MI_state_analysis(activity_matrix, states_matrix)
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

numNeurons = size(activity_matrix,2);
numStates = size(states_matrix,2);

MI = zeros(numStates,numNeurons);
p_values = zeros(numStates,numNeurons);
zMI = zeros(numStates,numNeurons);

for cell_i = 1:numNeurons
   for state_i = 1:numStates
   [MI(state_i,cell_i), p_values(state_i,cell_i), zMI(state_i,cell_i)] = compute_MI_significance(activity_matrix(:,cell_i), states_matrix(:,state_i));
   end
end

end

