function options = default_options_fit_vox_spectra(kernel,imgfilename)

if nargin == 1
    imgfilename = [];
end

options.kernel = kernel;


%% ILT options
options.ILT = default_ILT_options(kernel);


%% voxelwise integration options
%initialisation of the voxelwise clustering
%number of intergration regions (i.e. spectra ROIs)
options.nSROI = 4;

%default is to split the 
spectra_dim = length(options.ILT.Nk);
mink = options.ILT.mink;
maxk = options.ILT.maxk;

options.sROIbounds = cell(1,spectra_dim);

%calculate the number of sROI boundaries in each dimension 
SROI_dim = zeros(spectra_dim,1);
for i=1:spectra_dim-1
    %if the number of spectral dimensions doesn't divide the number
    %of sROI exactly, then take the ceiling of all dimensions ...
    SROI_dim(i) = ceil(options.nSROI / spectra_dim);
end
%.. except the last - take the floor
SROI_dim(spectra_dim) = floor(options.nSROI / spectra_dim);
   
for j=1:spectra_dim %loop through the kernel parameters
    %get the size of the sROIs for each parameter
    size_SROI(j) = (maxk(j) - mink(j))/SROI_dim(j);
end


for j=1:spectra_dim
    for i=1:SROI_dim(j)
        %set the bounds for each kernel dimension
        options.sROIbounds{j}(i,:) = ...
            [mink(j)+(i-1)*size_SROI(j) mink(j)+i*size_SROI(j)];
    end
end

%now make all combinations of the sROI boundaries
allgridcombs = ( combvec(options.sROIbounds{:}) )';

options.allgridcombs = allgridcombs;

%split the sROI boundaries - one cell for each kernel dim
for j=1:spectra_dim
    options.sROI{j} = options.allgridcombs(:,(2*j-1) : 2*j);
end









%% saving options
%flag
options.save = 1;

%directory name
options.dirname = ['vox_spectra' ...
    '_nSROI_' num2str(options.nSROI)  ...    
    '_alpha_' num2str(options.ILT.alpha)];    

%where to save the directory
options.save_path = [pwd '/'];
   
%basename for output filenames 
if ~isempty(imgfilename) %if input is filepath
    %get the image filename in a nice format
    if ~iscell(imgfilename)
        %remove any full path stuff
        [~,imgfilename] = fileparts(imgfilename);
        %remove nifti extensions
        options.scan_names{1} = remove_ext_from_nifti(imgfilename);
    else
        for i=1:length(imgfilename) %loop over all filenames
            %remove any full path stuff
            [~,imgfilename{i}] = fileparts(imgfilename{i});
            %remove nifti extensions
            options.scan_names{i} = remove_ext_from_nifti(imgfilename{i});
        end
    end
else %if input is array
    options.scan_names{1} = 'img';
end

    




