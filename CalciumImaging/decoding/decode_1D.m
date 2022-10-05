function [actual_bin, decoded_bin, mean_decoding_error, decoding_agreement] = decode_1D(binarized_data, behav_vec, bin_vec, cell_used, training_ts, decoding_ts)
%DECODE_1D Decode 1-dimensional variable
%   Takes in neuronal activity and behavior as inputs, and outputs decoding
%   error, agreement, etc
binSize = bin_vec(2) - bin_vec(1);
bin_centers_vec = bin_vec + binSize/2;
bin_centers_vec(end) = [];

numNeurons = size(binarized_data,2);

%% Extract tuning curves
for cell_i = 1:numNeurons
    [~, ~, occupancy_vector, prob_being_active(cell_i), tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), behav_vec, bin_vec, training_ts);
end

%% Measuring error
% Minimal a priori (use to remove experimental a priori)
occupancy_vector = occupancy_vector./occupancy_vector*(1/length(occupancy_vector));

[decoded_probabilities] = bayesian_decode1D(binarized_data, occupancy_vector, prob_being_active, tuning_curve_data, 1:numNeurons);

%% 1 - Spatial: Let us now estimate the mouse location using the maximum a posteriori (MAP) value
[max_decoded_prob, decoded_bin] = max(decoded_probabilities,[],1);
decoded_val = bin_centers_vec(decoded_bin);

% Before looking at the error rate, we must first bin the actual data using the same bin vector used by
% the decoder
actual_bin = nan*behav_vec;
actual_val = nan*behav_vec;
for bin_i = 1:length(bin_vec)-1
    position_idx = find(behav_vec>bin_vec(bin_i) & behav_vec < bin_vec(bin_i+1));
    actual_bin(position_idx) = bin_i;
    actual_val(position_idx) = bin_centers_vec(bin_i);
end

decoded_bin(~decoding_ts) = nan;
decoded_val(~decoding_ts) = nan;
decoded_probabilities(:,~decoding_ts) = nan;

actual_bin(~decoding_ts) = nan;
actual_val(~decoding_ts) = nan;
actual_bin = actual_bin';
actual_val =  actual_val';

decoding_agreement_vector = double(decoded_bin == actual_bin);
decoding_agreement_vector(isnan(decoded_bin)) = nan;
decoding_agreement_vector(isnan(actual_bin)) = nan;
decoding_agreement_vector(isnan(decoding_agreement_vector)) = [];
decoding_agreement = sum(decoding_agreement_vector)./length(decoding_agreement_vector);

decoding_error = actual_val - decoded_val;
mean_decoding_error = mean(abs(decoding_error), 'omitnan');

end

