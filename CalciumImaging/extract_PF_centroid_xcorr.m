for day=5:10
    dayPF{day-4} = [PV988{day} PV989{day} PV990{day} PV991{day}]; 
end

for cell_i=1:size(dayPF{1},2)
   ref = dayPF{1}(:,cell_i);
   for day_i = 1:length(dayPF)
      [coef(:,day_i,cell_i), lags] = xcorr(ref,dayPF{day_i}(:,cell_i),'coeff');
   end
end

binSize = 3;
lags = lags*binSize;

mean_coefs = mean(coef,3,'omitnan');