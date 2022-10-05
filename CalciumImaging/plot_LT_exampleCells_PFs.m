for day=5:10
    dayPF{day-4} = [PV986{day} PV988{day} PV989{day} PV990{day} PV991{day}]; 
end

[~, maxIdx] = max(dayPF{1});
[~, sortedIdx] = sort(maxIdx,'ascend');

figure
for day = 1:6
   data = dayPF{day};
   subplot(1,6,day)
   imagesc(data(:,sortedIdx)');
   daspect([1 1 1])
   colormap 'Viridis'
end