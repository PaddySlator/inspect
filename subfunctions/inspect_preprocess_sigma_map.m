function allsigma = inspect_preprocess_sigma_map(sigma,mask)

% preprocess sigma maps for inspect  
% 1. load sigma map files if necessary
% 2. concatenate multiple sigma maps into single array
% 
% sigma is either a single image volume, or multiple image volumes
% mask is output from inspect_processes_img


%check if sigma is a path to a nifti file/files
if ischar(sigma) || isstring(sigma) %single nifti file
    sigmafilename = char(sigma); %if string convert to character array
    %check if there are any spaces in the character array - if so, assume
    %that each space is separating a filename
    if sum(isspace(sigmafilename)) > 0 
        %split into cell array for loading later
        sigma=split(sigmafilename);
    else
        sigma = niftiread(sigmafilename); %load main image
        sigma = double(sigma);
    end
end
if iscell(sigma)
    if ischar(sigma{1}) || isstring(sigma{1}) %multiple nifti filenames in a cell      
        sigmafilename = cellfun(@char,sigma,'UniformOutput',false); %make sure it's a character
        %total number of filenames provided - i.e. number of images
        nimg = length(sigma);
        %store the images in a cell
        sigma = cell(nimg,1);
        for i=1:nimg
            sigma{i} = niftiread(sigmafilename{i});            
            sigma{i} = double(sigma{i});
        end
    else %multiple images in a cell
        sigmafilename=[];
    end
else
    sigmafilename=[];
end


%if the input SNR is a single image volume put it into a one element cell, then can
%just treat it the same as a multiple image stucture
if ~iscell(sigma)
    tempsigma = cell(1);

    tempsigma{1} = sigma;

    sigma =  tempsigma;
end


%%% concatenate all the SNR maps into a single image

%total number of images
nimg = length(sigma);

%for storing voxel positions - to map back into image space
sigmavox = cell(nimg,1);
voxind = cell(nimg,1);

allsigma = [];
imgind = [];

for i=1:nimg
    %put image into voxel form
    [sigmavox{i}, voxind{i}] = image_to_voxel(sigma{i},mask{i});
    allsigma = [allsigma; sigmavox{i}];
    %store the img index for these voxels
    imgind = [imgind; i*ones(size(sigmavox{i},1),1)];
end


