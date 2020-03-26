function mean_signal = calculate_mean_signal(dw_image,mask)
%calculate the mean signal across the voxels specified by a mask
%
%INPUTS
%dw_image - diffusion weighted image as a 4D array 
%mask - the mask to use to calculate the mean as a 3D array
%
%OUTPUTS
%mean_signal - the mean dw signal in the mask - has the same length as the number of diffusion volumes
%
%Author 
%Paddy Slator (p.slator@ucl.ac.uk)

image_size=size(dw_image(:,:,:,1));
n_voxels=nnz(mask(:,:,:,1));


signal=zeros(n_voxels,size(dw_image,4));

l=1;
for i=1:image_size(1)
    for j=1:image_size(2)
        for k=1:image_size(3)
            if mask(i,j,k)                
                signal(l,:)=dw_image(i,j,k,:);
                l=l+1;
            end
        end
    end
end




mean_signal=mean(signal,1);


end