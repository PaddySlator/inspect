function [voxel,voxel_index] = image_to_voxel(image,mask)
% convert from image coordinates to voxel coordinates

if nargin < 2
   mask = ones(size(image(:,:,:,1))); 
end

n_vol = size(image,4); 
n_voxels = nnz(mask);

voxel = zeros(n_voxels,n_vol);

%slow looped code
image_size=size(image);
l=1;
for i=1:image_size(1)
    for j=1:image_size(2)
        for k=1:image_size(3)            
            if mask(i,j,k)
                voxel(l,:)=squeeze(image(i,j,k,:));
                
                voxel_index(l,:)=[i j k];
                
                l=l+1;
            end
        end
    end
end


% vector=image(:);
% vector_mask=mask(:);
% 
% %remove the voxels not included in the mask
% vector=image(mask~=0);

end