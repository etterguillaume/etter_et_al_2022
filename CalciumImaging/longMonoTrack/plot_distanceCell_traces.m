function plot_distanceCell_traces(ms,behav,cellID)
%PLOT_TIMECELL_TRACES Summary of this function goes here
%   Detailed explanation goes here

%Params
lap_detect = 10;
spatialBinSize = 3;
max_cum_spatial_length = 536;

%% Extract variables from data structure
ca_time = ms.time/1000;
ca_trace = ms.RawTraces(:,cellID);
ca_trace=detrend(ca_trace);

behav_vec = behav.position(:,1);
behav_time = behav.time/1000;
 
optosignal = behav.optosignal;

%% Extract tone state
LMT_state = extractLMTState(optosignal);

[behav_time, IAbehav, ~] = unique(behav_time);
[ca_time, IAms, ~] = unique(ca_time);
ca_trace = ca_trace(IAms,:);
behav_vec = behav_vec(IAbehav,:);
LMT_state = LMT_state(IAbehav);

peakCa = max(ca_trace);

%% Interpolate behavior
interp_behav_vec = interpolate_behavior(behav_vec, behav_time, ca_time); % in the X dimension
interp_behav_vec(end) = interp_behav_vec(end-1);

%% Extract state
LMT_state = interp1(behav_time, LMT_state,ca_time,'nearest');

%% Extract traveled distance
travelled_dist = diff(interp_behav_vec);
travelled_dist(end+1)=0;
cum_distance(1) = 0;
cum_time(1) = 0;

traceNum = 0;
ctr=1;

tracePlot = [];

for frame_i = 2:length(LMT_state)
    if LMT_state(frame_i) == 1 && interp_behav_vec(frame_i-1) < lap_detect && interp_behav_vec(frame_i) >= lap_detect
        cum_distance(frame_i) = 0;
        traceNum = traceNum+1;
        ctr=1;
    end
    
    if traceNum>0 & cum_distance <= 536
        cum_distance(frame_i) = cum_distance(frame_i-1)+abs(travelled_dist(frame_i));
        tracePlot(traceNum,ctr) = ca_trace(frame_i);
        ctr=ctr+1;
    end
end

%tracePlot = tracePlot-min(tracePlot,[],2);
%tracePlot=tracePlot./max(tracePlot,[],2);

%traces2remove = [3 8 10 11];
%tracePlot(traces2remove,:) = [];
    
figure
for trace=1:traceNum-length(traces2remove)
   plot(timeScale,tracePlot(trace,:)+trace-1); hold on 
end
ax=gca;
ax.XLim = [0 20];

figure
imagesc(tracePlot)
ax2=gca;
ax.XLim = [0 20];

end

