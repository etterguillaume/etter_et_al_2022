dayPF={};
for day=1:10
    dayPF{day} = [PV986{day} PV988{day} PV990{day} PV991{day}]; 
end

figure
for day = 1:10
   data = dayPF{day};
   [~, maxIdx] = max(data);
   [~, sortedIdx] = sort(maxIdx,'descend');
   subplot(1,10,day)
   imagesc(data(:,sortedIdx)');
   daspect([1 1 1])
   colormap 'Viridis'
end

%% Sort to day 1
[~, maxIdx] = max(dayPF{1});
[~, sortedIdx] = sort(maxIdx,'descend');  
figure
for day = 1:10
   data = dayPF{day};
   subplot(1,10,day)
   imagesc(data(:,sortedIdx)');
   daspect([1 1 1])
   colormap 'Viridis'
end

%% Sort to day 10
[~, maxIdx] = max(dayPF{10});
[~, sortedIdx] = sort(maxIdx,'descend');      
figure
for day = 1:10
   data = dayPF{day};
   subplot(1,10,day)
   imagesc(data(:,sortedIdx)');
   daspect([1 1 1])
   colormap 'Viridis'
end