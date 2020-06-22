function [img,imgfilename,mask,nimg,allimg,imgind,voxind,nvox,nx,ny,nz] = inspect_preprocess_img(img,mask)

% preprocess the inspect input 
% 1. load image, mask and protocol ('gradechoinv') files if necessary
% 2. concatenate multiple images/masks
% 
% img is either a single image volume, or multiple image volumes


%check if img is a path to a nifti file/files
if ischar(img) || isstring(img) %single nifti file
    imgfilename = char(img); %if string convert to character
    img = niftiread(imgfilename); %load main image
elseif iscell(img)
    if ischar(img{1}) || isstring(img{1}) %multiple nifti filenames      
        imgfilename = cellfun(@char,img,'UniformOutput',false); %make sure it's a character
        %total number of filenames provided - i.e. number of images
        nimg = length(img);
        %store the images in a cell
        img = cell(nimg,1);
        for i=1:nimg
            img{i} = niftiread(imgfilename{i});
        end
    else %multiple images in a cell
        imgfilename=[];
    end
else
    imgfilename=[];
end

%check if mask is a path to a nifti file
if ischar(mask) || isstring(mask)
    maskfilename = mask;
    mask = niftiread(maskfilename); %load mask image
    
elseif iscell(mask)
    if ischar(mask{1}) || isstring(mask{1}) %multiple nifti filenames
        maskfilenames = mask;
        %store the images in a cell
        mask = cell(nimg,1);
        for i=1:nimg
            mask{i} = niftiread(maskfilenames{i});
        end
    end
else
    maskfilename=[];
end


%if the input is a single image volume put it into a one element cell, then can
%just treat it the same as a multiple image stucture
if ~iscell(img)
    tempimg = cell(1);
    tempmask = cell(1);

    tempimg{1} = img;
    tempmask{1} = mask;

    img =  tempimg;
    mask = tempmask;
end

%%% concatenate all the images into a single image

%total number of images
nimg = length(img);

%for storing voxel positions - to map back into image space
imgvox = cell(nimg,1);
voxind = cell(nimg,1);

allimg = [];
imgind = [];

for i=1:nimg
    %put image into voxel form
    [imgvox{i}, voxind{i}] = image_to_voxel(img{i},mask{i});
    allimg = [allimg; imgvox{i}];
    %store the img index for these voxels
    imgind = [imgind; i*ones(size(imgvox{i},1),1)];
end

%get the dimensions of the concatenated image
nvox = size(allimg,1);

%get the dimensions of each image
nx = cell(nimg,1);
ny = cell(nimg,1);
nz = cell(nimg,1);
for i=1:nimg
    %dimensions of each volume
    [nx{i},ny{i},nz{i},~] = size(img{i}(:,:,:,1));
end



