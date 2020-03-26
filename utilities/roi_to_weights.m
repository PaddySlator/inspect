function weights = roi_to_weights(roi,mask)
%convert weights (i.e. posterior probabilities) to a ML roi 


nclus = max(roi(:));

[nx,ny,nz] = size(roi);

if nargin < 2
   mask = ones([nx ny nz]); 
end

weights = zeros([nx ny nz nclus]);


for x=1:nx
    for y=1:ny
        for z=1:nz
            if mask(x,y,z)
                %find which cluster has the maximum weight for this voxel
                weights(x,y,z,roi(x,y,z)) = 1;                         
            end
        end
    end
end


end