function options = default_options_inspect_map(kernel,gradechoinv,imgfilename)

if nargin == 2
    imgfilename = [];
end

options.kernel = kernel;

%% algorithm options
%initialisation of the voxelwise clustering
options.init='meanspectrum';
%number of compartments
options.ncomp = 3;
%maximum iterations
options.maxiter = 5;
%tolerance for convergence
options.weightstol = 10^-3;
%turn updates on/off for debugging
options.onF = 1;
options.onweights = 1;
%choose whether to order the fitted components 
options.relabel = 1;
%estimate SNR by taking standard deviation of voxels
options.SNR = 'voxelwise';
%this is the default SNR value that will be used if cannot be estimated by
%taking standard deviation
options.fixedSNR = 50;
%not using a fixed sigma map so set to NaN
options.fixedsigma = NaN;
%choose whether to normalise the volume fractions
options.sumto1 = 1;


%% ILT options
options.ILT = default_ILT_options(kernel,gradechoinv);
%for the initial fit to the mean signal
options.ILT_mean = default_ILT_options(kernel,gradechoinv);



%% saving options
%flag
options.save = 1;

%directory name
options.dirname = ['inspect_map' ...
    '_ncomp_' num2str(options.ncomp)  ...
    '_maxiter_' num2str(options.maxiter) ...
    '_alpha_' num2str(options.ILT.alpha) ...
    '_init_' num2str(options.init)];    

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