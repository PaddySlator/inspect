function weighted_mean_signal = calculate_weighted_mean_signal(img,mask,weights)
%calculate the mean signal across the voxels specified by a mask, with the
%weight associated for each voxel
%
%INPUTS
%img -  image as a 4D array 
%mask - the mask to use to calculate the mean as a 3D array
%weights - an image with the same size as mask that specifies the weights
%to assign to each voxel in the mean signal
%
%OUTPUTS
%mean_signal - the mean dw signal in the mask - has the same length as the number of diffusion volumes
%
%Author 
%Paddy Slator (p.slator@ucl.ac.uk)

if isempty(mask)
    
end
    
   
image_size=size(img(:,:,:,1));
n_voxels=nnz(mask(:,:,:,1));


signal=zeros(n_voxels,size(img,4));




l=1;
for i=1:image_size(1)
    for j=1:image_size(2)
        for k=1:image_size(3)
            if mask(i,j,k)
                %weighted signal in this voxel
                signal(l,:)= (weights(i,j,k) * img(i,j,k,:))  ;
                l=l+1;
            end
        end
    end
end

%the sum of all the weights for normalising
sumweights = sum(sum(sum(weights.*mask)));
%final sum
weighted_mean_signal= sum(signal,1) / sumweights;
    
    
