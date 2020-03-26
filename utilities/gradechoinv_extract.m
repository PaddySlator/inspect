function gradechoinv_extract(img_filename,gradechoinv_filename,options)
%Extract volumes from diffusion-relaxometry image and grad echo inv file


%check if filename is an absolute path or just a filename
if isempty(fileparts(img_filename))
    %if just a filename then assume we are already in the right directory
    dwi_full_path=pwd;
    %split the filename into name and ext
    name=remove_ext_from_nifti(img_filename);
    %remove the extension from grad filename
    [~,grad_name]=fileparts(gradechoinv_filename);
else
    %absolute path given, so move to the directory
    cd(fileparts(img_filename))
    dwi_full_path=pwd;
    %get the index of the filename (part_name might not have the full extension)
    [~,part_name,~]=fileparts(img_filename);
    name_index=strfind(img_filename,part_name);
    %split the filename into name and ext
    name=remove_ext_from_nifti(img_filename(name_index:end));
    %remove the extension from grad filename
    [~,grad_name]=fileparts(gradechoinv_filename);
end


%load the nifti image
full_dwi=load_untouch_nii(img_filename);
%load the gradient file
grads=load(gradechoinv_filename);
%first 3 columns are b-vectors
bvecs=grads(:,1:3)';
%4th column is b-values
bvals=grads(:,4)';
%5th column is echo times
te = grads(:,5)';


%round the b-values to nearest integer 
bvals=round(bvals);
%find the index of b0 volumes
b0_index=find(bvals==0);
    
    

if ~isfield(options, 'shell')
    options.shell=0;
end
if ~isfield(options, 'abovecutoff')
    options.abovecutoff=0;
end
if ~isfield(options, 'belowcutoff')
    options.belowcutoff=0;
end
if ~isfield(options, 'volumes')
    options.volumes=0;
end


if options.shell %extract the given shells 
    %extract the index of the volumes to keep 
    index_to_keep{1}=b0_index;
    for i=1:length(options.shells)
        index_to_keep{1}=[index_to_keep{1} find(bvals == options.shells(i))];
    end
    %get the base filename for the nifti and grad files
    sub_dwi_name{1}=[name '_' replace_space_with_underscore(['b_' num2str(options.shells)])];
    sub_grad_name{1}=[grad_name '_' replace_space_with_underscore(['b_' num2str(options.shells)])];
else
    index_to_keep{1}=[];
    sub_dwi_name{1}=[];
end

if options.abovecutoff %extract all images above a given b-value   
    %extract the index of the volumes to keep 
    index_to_keep{2}=[b0_index find(bvals >= options.bmin)];
    %get the base filename for the nifti and grad files
    sub_dwi_name{2}=[name '_b_geq_' num2str(options.bmin)];
    sub_grad_name{2}=[grad_name '_b_geq_' num2str(options.bmin)];
else
    index_to_keep{2}=[];
    sub_dwi_name{2}=[];
end

if options.belowcutoff %extract all images below a given b-value
    %extract the index of the volumes to keep 
    index_to_keep{3}=find(bvals <= options.bmax);    
    %get the base filename for the nifti and grad files
    sub_dwi_name{3}=[name '_b_leq_' num2str(options.bmax)];
    sub_grad_name{3}=[grad_name '_b_leq_' num2str(options.bmax)];
else
    index_to_keep{3}=[];
    sub_dwi_name{3}=[];
end

if options.volumes %remove the chosen volumes
    index_to_keep{4} = options.index_to_keep;
    sub_dwi_name{4} = [name 'trunc'];
    sub_grad_name{4} = [grad_name 'trunc'];
else
    index_to_keep{4}=[];
    sub_dwi_name{4}=[];
end


for i=1:length(sub_dwi_name) %loop through the different subsamplings
    if ~isempty(sub_dwi_name{i})
        %extract the volumes from the full nifti file
        sub_dwi=full_dwi;
        %subsample image
        sub_dwi.img=full_dwi.img(:,:,:,index_to_keep{i});
        %change the dimensions in the header
        sub_dwi.hdr.dime.dim(2:5)=size(sub_dwi.img);
        %save as nifti file
        save_untouch_nii(sub_dwi,[sub_dwi_name{i} '.nii.gz'])
        
        %make the grad file
        %create txt file with write access (same directory as original nifti file)
        fid=fopen([sub_grad_name{i} '.txt'],'w');
        %get the rows of the grad table      
        grads(index_to_keep{i},:)
        sub_grad_file=grads(index_to_keep{i},:);
        %print to file
        fprintf(fid,'%4f %4f %4f %4f %4f %4f \n',sub_grad_file');
        %close the file
        fclose(fid);
                          

    end
end
    


    
    
    
    
    

end