canvas = zeros(size(cell_registered_struct.spatial_footprints_corrected{1},2),size(cell_registered_struct.spatial_footprints_corrected{1},3));

for cell_i = 1:size(cell_registered_struct.cell_to_index_map,1);
    day2use=1;
    while cell_registered_struct.cell_to_index_map(cell_i,day2use) == 0
        day2use=day2use+1;
    end
    canvas = canvas+permute(cell_registered_struct.spatial_footprints_corrected{day2use}(cell_registered_struct.cell_to_index_map(cell_i,day2use),:,:),[2 3 1]);
    
end

figure
pcolor(canvas)
shading interp
colormap Viridis
daspect([1 1 1]);
