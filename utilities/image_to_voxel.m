function [voxelimg,voxelimg_index] = image_to_voxel(image,mask)
% convert from image coordinates to voxel coordinates

if nargin < 2
   mask = ones(size(image(:,:,:,1))); 
end

n_vol = size(image,4); 
n_voxels = nnz(mask);


% tic;
% voxelimg = zeros(n_voxels,n_vol);
% voxelimg_index = zeros(n_voxels,3);
% 
% %slow looped code
% image_size=size(image);
% l=1;
% for k=1:image_size(3)
%     for j=1:image_size(2)
%         for i=1:image_size(1)            
%             if mask(i,j,k)
%                 voxelimg(l,:)=squeeze(image(i,j,k,:));
%                 
%                 voxelimg_index(l,:)=[i j k];
%                                 
%                 l=l+1;
%             end
%         end
%     end
% end
% toc;


tic;
%fast vectorised code
%get the indices of the non-zero mask voxels
[x,y,z] = ind2sub(size(mask),find(mask));
voxelimg_index = [x y z];

%extract these voxels from the image
%repeat the mask across all volumes
mask_allvol = repmat(mask,[1 1 1 n_vol]);
voxelimg = image(logical(mask_allvol));
voxelimg = reshape(voxelimg,[n_voxels n_vol]);
toc;




% vector=image(:);
% vector_mask=mask(:);
% 
% %remove the voxels not included in the mask
% vector=image(mask~=0);

end