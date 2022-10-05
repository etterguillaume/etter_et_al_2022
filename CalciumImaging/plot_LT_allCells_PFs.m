% for day=5:10
%     dayPF{day-4} = [PV988{day} PV989{day} PV990{day} PV991{day}];
% end


for cell_i = 1:size(dayPF{1},2)
    figure
    for day=1:6
        hold on
        if day == 3
            plot(dayPF{day}(:,cell_i)+day-1,'color', 'r')
        elseif day == 4
            plot(dayPF{day}(:,cell_i)+day-1,'color', 'b')
        else
            plot(dayPF{day}(:,cell_i)+day-1,'color', 'k')
        end
    end
    hold off
    drawnow
    pause
end