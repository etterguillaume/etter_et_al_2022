for day=5:10
    dayPF{day-4} = [PV988{day} PV989{day} PV990{day} PV991{day}]; 
end

for day_i = 1:length(dayPF)
    [maxVal, maxIdx(day_i,:)] = max(dayPF{day_i},[],1,'omitnan');
end

binSize = 3;
maxIdx = maxIdx*3;

centroidShift = maxIdx-maxIdx(1,:);

%% Compute absolute shift from day 1
absShift=abs(centroidShift);
mean(absShift,2);