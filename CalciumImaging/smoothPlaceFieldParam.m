function smoothed_map = smooth_map(input_map,smooth_factor)
   
if size(input_map,2) > 1 
    % Used to smooth place fields
    kernel_size = [3 3];
    occupancy_std = 2;
    [Xgrid,Ygrid]=meshgrid(-kernel_size(1)/2: kernel_size(1)/2, -kernel_size(2)/2:kernel_size(2)/2);
    Rgrid=sqrt((Xgrid.^2+Ygrid.^2));
    kernel = pdf('Normal', Rgrid, 0, occupancy_std);
    kernel = kernel./sum(sum(kernel));
    smoothed_map = conv2(input_map, kernel, 'same'); % smoothing 

else
   smoothed_map = smoothdata(input_map,'g',smooth_factor); % Default 5
end

end

