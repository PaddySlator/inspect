function options = default_options_inspect_seg(kernel,gradechoinv,imgfilename)

if nargin == 2
    imgfilename = [];
end

options.kernel = kernel;

%% EM algorithm options
%initialisation of the voxelwise clustering
options.init='kmeans';
%number of clusters
options.nclus = 4;
%maximum number of EM steps
options.maxiter = 10;
%tolerance for convergence
options.weightstol = 10^-3;



%% ILT (i.e. M-step) options
options.ILT = default_ILT_options(kernel,gradechoinv);


%% saving options
%flag
options.save = 1;

%directory name
options.dirname = ['inspect_seg' ...
    '_nclus_' num2str(options.nclus)  ...
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

    




