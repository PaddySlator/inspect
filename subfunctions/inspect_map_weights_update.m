function weights = inspect_map_weights_update(allimg,currentweights,Fcomp,Kalpha,gradechoinv,options)
% 
% Weight update step for inspect continuous mapping version.
% Given canonical spectra, estimates voxelwise spectral weightings.

% INPUTS:
%
% allimg - image data in voxelwise form, typically output from 
% inspect_preprocess_img.m
%
% currentweights - current values for the voxelwise weights
% to use as starting points for the optimisation. Pass "[]" to
% use random starting points.
%
% Fcomp - current canonical spectra  
%
% Kalpha - the dictionary of kernal values, including regularisation
% elements
%
% gradechoinv - array  of the MR acquisition parameters.  
% Format: [gx gy gz b TE TI TR] where gx gy gz are the gradient 
% directions, b is the b-value, TE is the echo time, TI is the inversion 
% time, and TR is the repetition time.  
%
% options - algorithm options
%
%
%
% Author: Paddy Slator, p.slator@ucl.ac.uk
%
%
% LICENSE
% <inspect toolbox for qMRI analysis>
% Copyright (C) <2020>  <Paddy J. Slator, p.slator@ucl.ac.uk>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.







%unpack some required parameters
Nk = options.ILT.Nk;
reg = options.ILT.reg;
ncomp = options.ncomp;
nvox = size(allimg,1);

%get 1/SNR for all voxels
sigvec = estimate_sd(allimg,gradechoinv);  

%get starting point for weights (use current weight if available)
if isempty(currentweights)
    z = rand(nvox,ncomp);
    %normalise
    z = z./repmat(sum(z,2),[1 ncomp]);
else
    z = currentweights;
end
    
%set up fmincon options
fmincon_options = optimoptions(@fmincon,'Display','off');

%set up lower bounds and upper bounds - make sure fitted params are probabilities
lb = zeros(1,ncomp);
ub = ones(1,ncomp);



%preallocate weights
weights = zeros(nvox,ncomp);

parfor i=1:nvox
    ...for i=1:nvox    
        %get the fmincon starting point 
        z0 = z(i,:);
        %double check that they are normalised to 1 - and hence satisfy the
        %constraint
        z0 = z0./sum(z0);
        
        %get the signal for this voxel
        S = allimg(i,:)';
        
        if reg
            %augment the signal
            S = [S; zeros(prod(Nk), 1)];
        end
                             
        %make a function for optimisation (minimize the negative logli) which takes only weights as input
        zoptfun = @(z) -calculate_map_logli(S,z,Fcomp,sigvec(i),Kalpha);
        
        %optimize for weights, including probcon function (parameters sum to 1)
        weights(i,:) = fmincon(zoptfun,z0,[],[],[],[],lb,ub,@probcon,fmincon_options);
        
        nupdate=1;%display an update for every nth voxel
        if mod(i,nupdate) == 0
            disp(['voxel ' num2str(i) ' of ' num2str(nvox)])
        end               
end
    

    
    

                  
        
end