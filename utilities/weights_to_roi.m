function MLroi = weights_to_roi(weights,mask)
%convert weights (i.e. posterior probabilities) to a ML roi 

[nx,ny,nz,nclus] = size(weights);

if nargin < 2
   mask = ones([nx ny nz]); 
end

MLroi = zeros([nx ny nz]);

for x=1:nx
    for y=1:ny
        for z=1:nz
            if mask(x,y,z)
                %find which cluster has the maximum weight for this voxel
                [~,maxind] = max(weights(x,y,z,:));
                
                MLroi(x,y,z) = maxind;
            end
        end
    end
end






end