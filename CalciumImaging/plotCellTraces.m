function plotCellTraces(traces)
%MSPLOTCELLTRACEAMONGOTHERS Summary of this function goes here
%   Detailed explanation goes here

gain = 5;

time = (0:length(traces))/30;
time(end) = [];

figure
%for trace_i = 1:size(traces,2)
for trace_i = 1:100
        plot(time,traces(:,trace_i)*gain+2*(trace_i-1),'color',[0.8 0.8 0.8]); hold on;
end
ax=gca;
ax.XLim = [0 60];


end