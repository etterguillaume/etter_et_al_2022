%% Parameters
bootstrapSamples = 30;
binSize = 3;
sampling_frequency = 30; % This data set has been sampled at 30 images per second
z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise
min_speed_threshold = 5; % 2 cm.s-1
OpticSigThreshold = 0.2;

%% Create bins
X_bin_vector = 0:binSize:45+binSize;
X_bin_centers_vector = X_bin_vector + binSize/2;
X_bin_centers_vector(end) = [];

Y_bin_vector = 0:binSize:45+binSize;
Y_bin_centers_vector = Y_bin_vector + binSize/2;
Y_bin_centers_vector(end) = [];

%% Indicate working folder here
workingFolder = pwd;
split_path = strsplit(workingFolder,filesep);
mouse_number = split_path{end-1};

%% Analyze folder structure
session_dates = [];
session_folders = dir(workingFolder);
folder2process = [];
for file_i = 1:length(session_folders)
    split_session = strsplit(session_folders(file_i).name,'_');
    if strcmp(split_session{1},mouse_number) 
        session_dates(file_i) = str2num(split_session{end}); % For potential chronological sorting later
        folder2process(file_i) = 1;
    else
        session_dates(file_i) = 0;
        folder2process(file_i) = 0;
    end
    if strcmp(split_session{1},'cellRegistered') 
        load([workingFolder filesep session_folders(file_i).name]);
        cell_to_index_map = cell_registered_struct.cell_to_index_map; % Load the identity matrix
    end
end

[sorted_dates, sorted_idx] = sort(session_dates);

ca_data = {};
ca_time = {};
behav_vec = {};
behav_time = {};
optoSignal = {};

sesh = 1;
for folder_i = 1:length(sorted_dates)
    currentFolder = sorted_idx(folder_i);
    if folder2process(currentFolder)
        load([workingFolder filesep session_folders(currentFolder).name '/ms.mat'])
        ca_data{sesh} = ms.RawTraces;
		ca_time{sesh} = ms.time/1000;
		load([workingFolder filesep session_folders(currentFolder).name '/behav.mat'])
        behav_vec{sesh} = behav.position;
		behav_time{sesh} = behav.time/1000;

		if isfield(behav,'optosignal')
			optoSignal{sesh} = behav.optosignal;
		else
			optoSignal{sesh} = [];
		end
		sesh = sesh +1;
    end
end


%% Pre-allocate memory and create labels
withinDaysCorr = zeros(30,9);
withinDaysLabels = {'Within day 1', 'Within day 2', 'Within day 3 baseline', 'Between day 3 baseline & scrambled stims', 'Within day 3 scrambled stims', 'Within day 4 baseline', 'Between day 4 baseline & 8 Hz stims', 'Within day 4 8Hz stims', 'Within day 5'};

%% Preprocess the data
interp_behav_vec = {};
binarizedData = {};
running_ts = {};
stimulatedEpochs = {};
Fs = sampling_frequency;
%[bFilt,aFilt] = butter(2,  2/(30/2), 'low');

StimDuration = 5; % in s
StimDurationFrames = round(StimDuration*Fs);

for session_i = 1:length(ca_data)
	% Only keep unique time points
	[behav_time{session_i}, IAbehav, ~] = unique(behav_time{session_i});
	[ca_time{session_i}, IAms, ~] = unique(ca_time{session_i});
	ca_data{session_i} = ca_data{session_i}(IAms,:);
	behav_vec{session_i} = behav_vec{session_i}(IAbehav,:);
	numNeurons{session_i} = size(ca_data{session_i},2);
	
	%% Interpolate behavior
	interp_behav_vec{session_i}(:,1) = interpolate_behavior(behav_vec{session_i}(:,1), behav_time{session_i}, ca_time{session_i}); % in the X dimension
	interp_behav_vec{session_i}(:,2) = interpolate_behavior(behav_vec{session_i}(:,2), behav_time{session_i}, ca_time{session_i}); % in the Y dimension
	interp_behav_vec{session_i}(end,:) = interp_behav_vec{session_i}(end-1,:);

	%% Extract velocity
	velocity{session_i} = extract_velocity(interp_behav_vec{session_i}, ca_time{session_i});

	%% Binarize
	binarizedData{session_i} = 0*ca_data{session_i};
	for cell_i = 1:numNeurons{session_i}
	    binarizedData{session_i}(:,cell_i) = extract_binary(ca_data{session_i}(:,cell_i), sampling_frequency, z_threshold);
	end
	running_ts{session_i} = velocity{session_i} > min_speed_threshold;

	%% Process opto signal
	if ~isempty(optoSignal{session_i})
        optoSignal{session_i} = optoSignal{session_i}(IAbehav);
		interpolatedOptoSignal = interp1(behav_time{session_i},optoSignal{session_i},ca_time{session_i});
		smoothedOptoSignal = smooth(interpolatedOptoSignal, 'sgolay');
		normalizedOptoSignal = smoothedOptoSignal-min(smoothedOptoSignal);
		normalizedOptoSignal = normalizedOptoSignal./max(normalizedOptoSignal);

		stimulatedEpochs{session_i} = zeros(size(normalizedOptoSignal));
		stimulatedEpochs{session_i}(normalizedOptoSignal>OpticSigThreshold) = 1;
	end
end


%% Start the bootstrapping procedure
for bootstrap_i = 1:bootstrapSamples
	disp(['Bootstrapping:' num2str(round(bootstrap_i/bootstrapSamples*100)) '%'])
	%% Day 1-A,B
    sample_A = running_ts{1};
    sample_B = sample_A;
    for i = 1:length(sample_A)
		if sample_A(i) == 1 && rand < 0.5
		 	sample_A(i) = 0;
		else
	    	sample_B(i) = 0;
		end
    end

	for cell_i = 1:numNeurons{1}
	    [~, ~, ~, ~, PF] = extract_2D_information(binarizedData{1}(:,cell_i), interp_behav_vec{1}, X_bin_vector, Y_bin_vector, sample_A);
	    PF_D1A(:,:,cell_i) = smoothPlaceField(PF);
	    [~, ~, ~, ~, PF] = extract_2D_information(binarizedData{1}(:,cell_i), interp_behav_vec{1}, X_bin_vector, Y_bin_vector, sample_B);
	    PF_D1B(:,:,cell_i) = smoothPlaceField(PF);
	end

	%% Day 2-A,B
    sample_A = running_ts{2};
    sample_B = sample_A;
    for i = 1:length(sample_A)
		if sample_A(i) == 1 && rand < 0.5
		 	sample_A(i) = 0;
		else
	    	sample_B(i) = 0;
		end
    end

	for cell_i = 1:numNeurons{2}
	    [~, ~, ~, ~, PF] = extract_2D_information(binarizedData{2}(:,cell_i), interp_behav_vec{2}, X_bin_vector, Y_bin_vector, sample_A);
	    PF_D2A(:,:,cell_i) = smoothPlaceField(PF);
	    [~, ~, ~, ~, PF] = extract_2D_information(binarizedData{2}(:,cell_i), interp_behav_vec{2}, X_bin_vector, Y_bin_vector, sample_B);
	    PF_D2B(:,:,cell_i) = smoothPlaceField(PF);
	end

	%% Day 5-A,B
    sample_A = running_ts{5};
    sample_B = sample_A;
    for i = 1:length(sample_A)
		if sample_A(i) == 1 && rand < 0.5
		 	sample_A(i) = 0;
		else
	    	sample_B(i) = 0;
		end
    end

	for cell_i = 1:numNeurons{5}
	    [~, ~, ~, ~, PF] = extract_2D_information(binarizedData{5}(:,cell_i), interp_behav_vec{5}, X_bin_vector, Y_bin_vector, sample_A);
	    PF_D5A(:,:,cell_i) = smoothPlaceField(PF);
	    [~, ~, ~, ~, PF] = extract_2D_information(binarizedData{5}(:,cell_i), interp_behav_vec{5}, X_bin_vector, Y_bin_vector, sample_B);
	    PF_D5B(:,:,cell_i) = smoothPlaceField(PF);
	end

	%% Day 3_baseline-A,B
    sample_A = running_ts{3} & ~stimulatedEpochs{3}; %% Non-stimulated epochs
    sample_B = sample_A;
    for i = 1:length(sample_A)
		if sample_A(i) == 1 && rand < 0.5
		 	sample_A(i) = 0;
		else
	    	sample_B(i) = 0;
		end
    end

	for cell_i = 1:numNeurons{3}
	    [~, ~, ~, ~, PF] = extract_2D_information(binarizedData{3}(:,cell_i), interp_behav_vec{3}, X_bin_vector, Y_bin_vector, sample_A);
	    PF_D3_baseline_A(:,:,cell_i) = smoothPlaceField(PF);
	    [~, ~, ~, ~, PF] = extract_2D_information(binarizedData{3}(:,cell_i), interp_behav_vec{3}, X_bin_vector, Y_bin_vector, sample_B);
	    PF_D3_baseline_B(:,:,cell_i) = smoothPlaceField(PF);
	end

	%% Day 3_scrambled-A,B
    sample_A = running_ts{3} & stimulatedEpochs{3}; %% Stimulated epochs
    sample_B = sample_A;
    for i = 1:length(sample_A)
		if sample_A(i) == 1 && rand < 0.5
		 	sample_A(i) = 0;
		else
	    	sample_B(i) = 0;
		end
    end

	for cell_i = 1:numNeurons{3}
	    [~, ~, ~, ~, PF] = extract_2D_information(binarizedData{3}(:,cell_i), interp_behav_vec{3}, X_bin_vector, Y_bin_vector, sample_A);
	    PF_D3_scrambled_A(:,:,cell_i) = smoothPlaceField(PF);
	    [~, ~, ~, ~, PF] = extract_2D_information(binarizedData{3}(:,cell_i), interp_behav_vec{3}, X_bin_vector, Y_bin_vector, sample_B);
	    PF_D3_scrambled_B(:,:,cell_i) = smoothPlaceField(PF);
	end

	%% Day 4_baseline-A,B
    sample_A = running_ts{4} & ~stimulatedEpochs{4}; %% Non-stimulated epochs
    sample_B = sample_A;
    for i = 1:length(sample_A)
		if sample_A(i) == 1 && rand < 0.5
		 	sample_A(i) = 0;
		else
	    	sample_B(i) = 0;
		end
    end

	for cell_i = 1:numNeurons{4}
	    [~, ~, ~, ~, PF] = extract_2D_information(binarizedData{4}(:,cell_i), interp_behav_vec{4}, X_bin_vector, Y_bin_vector, sample_A);
	    PF_D4_baseline_A(:,:,cell_i) = smoothPlaceField(PF);
	    [~, ~, ~, ~, PF] = extract_2D_information(binarizedData{4}(:,cell_i), interp_behav_vec{4}, X_bin_vector, Y_bin_vector, sample_B);
	    PF_D4_baseline_B(:,:,cell_i) = smoothPlaceField(PF);
	end

	%% Day 4_8Hz-A,B
    sample_A = running_ts{4} & stimulatedEpochs{4}; %% Stimulated epochs
    sample_B = sample_A;
    for i = 1:length(sample_A)
		if sample_A(i) == 1 && rand < 0.5
		 	sample_A(i) = 0;
		else
	    	sample_B(i) = 0;
		end
    end

	for cell_i = 1:numNeurons{4}
	    [~, ~, ~, ~, PF] = extract_2D_information(binarizedData{4}(:,cell_i), interp_behav_vec{4}, X_bin_vector, Y_bin_vector, sample_A);
	    PF_D4_8Hz_A(:,:,cell_i) = smoothPlaceField(PF);
	    [~, ~, ~, ~, PF] = extract_2D_information(binarizedData{4}(:,cell_i), interp_behav_vec{4}, X_bin_vector, Y_bin_vector, sample_B);
	    PF_D4_8Hz_B(:,:,cell_i) = smoothPlaceField(PF);
    end

        %% TEMP CODE (for plotting)
    subplot(2,2,1)
    imagesc(PF_D3_baseline_A(:,:,2));
    ax=gca;
    ax.CLim=[0 0.1];
    daspect([1 1 1])
    title(['Mouse: ' mouse_number ', cell: 2, Day 3, Baseline A'])
    colormap Viridis
    colorbar
    subplot(2,2,2)
    imagesc(PF_D3_scrambled_A(:,:,2));
    ax=gca;
    ax.CLim=[0 0.1];
    daspect([1 1 1])
    title(['Mouse: ' mouse_number ', cell: 2, Day 3, Scrambled A'])
    colormap Viridis
    colorbar
    subplot(2,2,3)
    imagesc(PF_D3_baseline_B(:,:,2));
    ax=gca;
    ax.CLim=[0 0.1];
    daspect([1 1 1])
    title(['Mouse: ' mouse_number ', cell: 2, Day 3, Baseline B'])
    colormap Viridis
    colorbar
    subplot(2,2,4)
    imagesc(PF_D3_scrambled_B(:,:,2));
    ax=gca;
    ax.CLim=[0 0.1];
    daspect([1 1 1])
    title(['Mouse: ' mouse_number ', cell: 2, Day 3, Scrambled B'])
    colormap Viridis
    colorbar
    pause
    
%% Compute within session correlations
% Day 1
for cell_i = 1:numNeurons{1}
	vec_A = PF_D1A(:,:,cell_i);
	vec_B = PF_D1B(:,:,cell_i);
	vec_A = vec_A(:);
	vec_B = vec_B(:);
	temp_corr(cell_i) = corr(vec_A,vec_B);
end

withinDaysCorr(bootstrap_i,1) = mean(temp_corr,'omitnan');

% Day 2
for cell_i = 1:numNeurons{2}
	vec_A = PF_D2A(:,:,cell_i);
	vec_B = PF_D2B(:,:,cell_i);
	vec_A = vec_A(:);
	vec_B = vec_B(:);
	temp_corr(cell_i) = corr(vec_A,vec_B);
end

withinDaysCorr(bootstrap_i,2) = mean(temp_corr,'omitnan');

% Day 3 - baseline
for cell_i = 1:numNeurons{3}
	vec_A = PF_D3_baseline_A(:,:,cell_i);
	vec_B = PF_D3_baseline_B(:,:,cell_i);
	vec_A = vec_A(:);
	vec_B = vec_B(:);
	temp_corr(cell_i) = corr(vec_A,vec_B);
end

withinDaysCorr(bootstrap_i,3) = mean(temp_corr,'omitnan');

% Day 3 - baseline vs scrambled
for cell_i = 1:numNeurons{3}
	vec_A = PF_D3_baseline_A(:,:,cell_i);
	vec_B = PF_D3_scrambled_B(:,:,cell_i);
	vec_A = vec_A(:);
	vec_B = vec_B(:);
	temp_corr(cell_i) = corr(vec_A,vec_B);
end

withinDaysCorr(bootstrap_i,4) = mean(temp_corr,'omitnan');

% Day 3 - scrambled
for cell_i = 1:numNeurons{3}
	vec_A = PF_D3_scrambled_A(:,:,cell_i);
	vec_B = PF_D3_scrambled_B(:,:,cell_i);
	vec_A = vec_A(:);
	vec_B = vec_B(:);
	temp_corr(cell_i) = corr(vec_A,vec_B);
end

withinDaysCorr(bootstrap_i,5) = mean(temp_corr,'omitnan');

% Day 4 - baseline
for cell_i = 1:numNeurons{4}
	vec_A = PF_D4_baseline_A(:,:,cell_i);
	vec_B = PF_D4_baseline_B(:,:,cell_i);
	vec_A = vec_A(:);
	vec_B = vec_B(:);
	temp_corr(cell_i) = corr(vec_A,vec_B);
end

withinDaysCorr(bootstrap_i,6) = mean(temp_corr,'omitnan');

% Day 4 - baseline vs 8Hz
for cell_i = 1:numNeurons{4}
	vec_A = PF_D4_baseline_A(:,:,cell_i);
	vec_B = PF_D4_8Hz_B(:,:,cell_i);
	vec_A = vec_A(:);
	vec_B = vec_B(:);
	temp_corr(cell_i) = corr(vec_A,vec_B);
end

withinDaysCorr(bootstrap_i,7) = mean(temp_corr,'omitnan');

% Day 4 - 8Hz
for cell_i = 1:numNeurons{4}
	vec_A = PF_D4_8Hz_A(:,:,cell_i);
	vec_B = PF_D4_8Hz_B(:,:,cell_i);
	vec_A = vec_A(:);
	vec_B = vec_B(:);
	temp_corr(cell_i) = corr(vec_A,vec_B);
end

withinDaysCorr(bootstrap_i,8) = mean(temp_corr,'omitnan');

% Day 5
for cell_i = 1:numNeurons{5}
	vec_A = PF_D5A(:,:,cell_i);
	vec_B = PF_D5B(:,:,cell_i);
	vec_A = vec_A(:);
	vec_B = vec_B(:);
	temp_corr(cell_i) = corr(vec_A,vec_B);
end
withinDaysCorr(bootstrap_i,9) = mean(temp_corr,'omitnan');

end

save([workingFolder filesep 'withinDaysPFCorr.mat'],'withinDaysCorr', 'withinDaysLabels', 'mouse_number')




















