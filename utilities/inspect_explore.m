function output = inspect_explore(img,gradechoinv,mask,kernel,options)

%test different inspect models on one/a few scans 
%INPUTS
% gradechoinv - either a single gradechoinv file  



%TO DO - MAKE THIS WORK WITH ARRAY JOBS!


%%  single fits on each scan

%make sure that the relevant inputs are strings

%estimate voxelwise spectra and volume fraction estimation
if nargin < 5 %no options specified
    output.vox = inspect_vox(img,gradechoinv,mask,kernel);
else
    output.vox = inspect_vox(img,gradechoinv,mask,kernel,options.vox);
end

%fit continuous inspect
if nargin < 5
    output.map = inspect_map(img,gradechoinv,mask,kernel);
else
    output.map = inspect_map(img,gradechoinv,mask,kernel,options.map);
end


%fit clustering inspect
if nargin < 5
    output.seg = inspect_seg(img,gradechoinv,mask,kernel);
else
    output.seg = inspect_seg(img,gradechoinv,mask,kernel,options.seg);
end
    



end