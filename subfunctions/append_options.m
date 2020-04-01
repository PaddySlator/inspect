function appended_options = append_options(options,default_options)  

%get a list of all the legitimate options
legit_ILT_options = fields(default_options.ILT);
legit_options = fields(default_options);

%user options for the ILT
if isfield(options,'ILT')
    %get the options which have been user-specified
    user_ILT_options = fields(options.ILT);
    
    for i=1:length(user_ILT_options) %loop through replacing default values
        if any(strcmp(legit_ILT_options,user_ILT_options{i})) %check if this is a legitimate ILT option
            default_options.ILT.(user_ILT_options{i}) = options.ILT.(user_ILT_options{i});
            
            disp(['Using user defined ILT option: ' ...
                user_ILT_options{i} ' = ' num2str(options.ILT.(user_ILT_options{i}))])
        end
    end
end

%all other user options
user_options = fields(options);
user_options(strcmp(user_options,'ILT'))=[]; %already done these

for i=1:length(user_options) %loop through replacing default values
    if any(strcmp(legit_options,user_options{i})) %check if this is a legitimate option
        %if class(options.(user_options{i})) == class(default_options.(user_options{i})) %check if the user supplied value is in the right format
        default_options.(user_options{i}) =  options.(user_options{i});
        
        if isnumeric(options.(user_options{i}))
            user_option_str = num2str(options.(user_options{i}));
        else
            user_option_str = options.(user_options{i});
        end
        
        disp(['Using user defined option: ' ...
            user_options{i} ' = ' user_option_str ])
            
            
        %else
        %    disp(['User supplied value for ' user_options{i} ' not in correct format.'...
        %        ' It should be a ' class(default_options.(user_options{i})) '.' ...
        %        ' Using default value.'])
        %end
    end
end

%assign the updated default options
appended_options = default_options;