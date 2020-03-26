function image = voxel_to_image(voxel_vector,voxel_index,image_size)

image=zeros([image_size size(voxel_vector,2)]);

for i=1:size(voxel_vector,1)          
    image(voxel_index(i,1),voxel_index(i,2),voxel_index(i,3),:)=voxel_vector(i,:);    
end          


end