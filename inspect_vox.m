function [output,outputsummary] = inspect_vox(img,gradechoinv,mask,kernel,options)

% Fit spectrum voxelwise and calculate volume fraction maps 
% by spectral integration in user defined regions.
%
% Reference (this is the earliest one I have found):
% English,A.E. et al. Quantitative Two-Dimensional time 
% Correlation Relaxometry.  MRM 22(2), 425?434 (dec 1991).  
% https://doi.org/10.1002/mrm.1910220250
%
%
% REQUIRED INPUTS:
%
% img - either a single image, a cell of multiple images, or file path to a
% single nifti image. Each image has dimension [dimx dimy dimz encodings]. 
% If using multiple images these must all have the same MR encodings (i.e. gradechoinv file). 
% For multiple images the model is fit to all images at once.
%
% gradechoinv - array, or path to a .txt file, of the MR acquisition 
% parameters.  Format: [gx gy gz b TE TI TR] where gx gy gz are the gradient 
% directions, b is the b-value, TE is the echo time, TI is the inversion 
% time, and TR is the repetition time.
%
% mask - either a single binary image a cell of multiple binary images,
% or file path to a single nifti binary mask.
% Will be used to choose the voxels in which to fit the model. The mask
% must have dimension [dimx dimy dimz] - i.e. matching the first three
% dimensions of the corresponding image.
%
% kernel - string specifying the choice of kernel, which relates
% the MR acquisition parameters to the MR signal. For example:
% - 'DT2' for diffusion-T2, for which the kernel equation 
%   is exp(-b*D)*exp(-TE/t2). 
% - 'T1' for T1 inversion recovery, for which the kernel equation 
%   is abs(1-2*exp(-TI/t1) + exp(-TR/t1)).
%
% OPTIONAL INPUTS:
%
% options - algorithm options
%
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



%% preprocess the image/images

[~,imgfilename,mask,nimg,allimg,imgind,voxind,nvox,nx,ny,nz] = inspect_preprocess_img(img,mask);


%% extract the MR acqusition parameters 

if ischar(gradechoinv)%check if gradechoinv is a path to a file
    gradechoinvfilename = gradechoinv;
    gradechoinv = importdata(gradechoinvfilename);
end


%% unpack algorithm options
 
%get the default options
default_options = default_options_inspect_vox(kernel,gradechoinv,imgfilename); 

if nargin < 5 %if no user defined options    
    options=default_options;
    disp('Using default options.')                
elseif nargin == 5 %amend any user defined options
    options = append_options(options,default_options);
end

%print the options that are going to be used
options
options.ILT

if options.save
    %print save directory, and create it if it doesn't exist
    if exist([options.save_path options.dirname], 'dir')
        disp(['Results will be saved at: ' options.save_path options.dirname])
    else
        mkdir([options.save_path options.dirname])
        disp(['Created directory for saving results at: ' options.save_path options.dirname])
    end
    disp(['Output filenames will end in '  strjoin(options.scan_names,'_')])
end




%% get the kernel dictionary values by doing a dummy fit to the first voxel

%this turns off the actual ILT calculation
options.ILT.onILT = 0; 
%do the dummy ILT calculation
ILT_output_test = ILT(allimg(1,:)', gradechoinv, options.ILT);
%store the values
options.ILT.K = ILT_output_test.K;
options.ILT.grid = ILT_output_test.grid;
if options.ILT.reg
    options.ILT.Kalpha = ILT_output_test.Kalpha;
end
%normal fitting for all subsequent calls
options.ILT.onILT=1;


%% unpack a few well-used options


nSROI = size(options.sROI{1},1);


%%

%preallocate voxel-space spectral volume fraction image
spectvf = zeros([nvox nSROI]);

spectra = cell(nvox,1);

for i=1:nvox
    %calculate the spectra (takes up a lot of memory so replace each time)    
    voxILT = ILT(allimg(i,:),gradechoinv,options.ILT);
            
    %calculate volume fraction map
    spectvf(i,:) = integrate_vox_spectrum(voxILT.F,voxILT.grid,options.sROI);    
    
    spectra{i}=voxILT;
    
    disp(['completed voxel ' num2str(i) ' of ' num2str(nvox)])
end
   
for i=1:nimg
    vfimg{i} = voxel_to_image(spectvf(imgind == i,:),...
        voxind{i},...
        [nx{i} ny{i} nz{i}]);
end


%save the volume fraction image as nifti
if isfield(options,'save')
    if options.save
        for i=1:nimg %separate nifti for each image        
            niftiwrite(vfimg{i},[options.save_path options.dirname '/Vfvoxelwise_'  options.scan_names{i} '.nii.gz'])
        end
        output.fullsavepath = [options.save_path options.dirname];
    end
end



output.spectvf = spectvf;
output.vfimg = vfimg;

%save the ILT output for one voxel
output.voxILT = voxILT;

output.spectra = spectra;

output.kernel = kernel;
output.options = options;

