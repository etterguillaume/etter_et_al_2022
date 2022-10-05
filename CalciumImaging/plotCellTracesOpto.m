function plotCellTracesOpto(ms,neurons2plot, time_segment,behav)
%MSPLOTCELLTRACEAMONGOTHERS Summary of this function goes here
%   Detailed explanation goes here

gain = 5;

if isempty(time_segment)
    time_segment = [ms.time(1)/1000 ms.time(end)/1000];
end

if isempty(neurons2plot)
    neurons2plot = 1:size(ms.RawTraces,2)
else
    neurons2plot = neurons2plot(1):neurons2plot(end)
end


[unique_behavtime, IAbehav, ICbehav]=unique(behav.time);
[unique_mstime, IAms, ICms]=unique(ms.time);

unique_mstime = unique_mstime/1000;

optoSignal = behav.optosignal;
InterpolatedOptoSignal = interp1(behav.time(IAbehav),optoSignal(IAbehav),ms.time(IAms));


figure
subplot(2,1,1)
for trace_i = 1:length(neurons2plot)
        plot(unique_mstime,ms.RawTraces(:,neurons2plot(trace_i))*gain+2*(trace_i-1),'color',[0.8 0.8 0.8]); hold on;
end
ax1 = gca;

subplot(2,1,2)
plot(unique_mstime, InterpolatedOptoSignal);
ax2 = gca;

linkaxes([ax1 ax2], 'x')
ax1.XLim = time_segment;

%ax.YLim(1) = 0;
%ax.YLim(2) = max(ms.RawTraces(:,ms.numNeurons)*1+2*(ms.numNeurons-1));
end