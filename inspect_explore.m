function output = inspect_explore(img,gradechoinv,mask,kernel,ncomp)

%test different inspect models on one/a few scans 
%INPUTS
% gradechoinv - either a single gradechoinv file  




%TO DO - MAKE THIS WORK WITH ARRAY JOBS!


%%  single fits on each scan

%make sure that the relevant inputs are strings



%fit continuous inspect
if nargin < 5    
    output.map = inspect_map(img,gradechoinv,mask,kernel);
else
    options.map.ncomp = ncomp;
    output.map = inspect_map(img,gradechoinv,mask,kernel,options.map);
end


%fit clustering inspect
if nargin < 5
    output.seg = inspect_seg(img,gradechoinv,mask,kernel);
else
    options.seg.nclus = ncomp;
    output.seg = inspect_seg(img,gradechoinv,mask,kernel,options.seg);
end
    
%estimate voxelwise spectra and volume fraction estimation
%just use default options
output.vox = inspect_vox(img,gradechoinv,mask,kernel);






end