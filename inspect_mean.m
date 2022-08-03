function output = inspect_mean(img,gradechoinv,mask,kernel,options)



%% preprocess the image/images

[~,imgfilename,mask,nimg,allimg,imgind,voxind,nvox,nx,ny,nz] = inspect_preprocess_img(img,mask);


%% extract the MR acqusition parameters 

if ischar(gradechoinv)%check if gradechoinv is a path to a file
    gradechoinvfilename = gradechoinv;
    gradechoinv = importdata(gradechoinvfilename);
end



%% unpack algorithm options

%get the default options
default_options=default_options_inspect_mean(kernel,gradechoinv);

if nargin < 5 %if no user defined options    
    options=default_options;
    disp('Using default options.')                
elseif nargin == 5 %amend any user defined options
    options = append_options(options,default_options);
end

%%

%calculate the mean signal
meansig = mean(allimg)';

%calculate the ILT
output = ILT(meansig,gradechoinv,options.ILT);
%output the options too
output.options = options;
%can save the whole output as it shouldn't be too big
outputsummary = output;

%%
if isfield(options,'save')
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
end

%%
if isfield(options,'save')
    if options.save
        try
            save([options.save_path options.dirname '/outputsummary_' strjoin(options.scan_names,'_')],'outputsummary');
        catch %if filename too long
            try
                save([options.save_path options.dirname '/outputsummary_' options.scan_names{1} '_and_more'],'outputsummary');
                disp('filename was too long! check outputsummary.options.scan_names for scans')
            catch %if the file is bigger than 2GB - not ideal but still save the output
                save([options.save_path options.dirname '/outputsummary_' strjoin(options.scan_names,'_')],'outputsummary','-v7.3');
            end
        end
    end
end


