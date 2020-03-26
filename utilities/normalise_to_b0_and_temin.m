function normimg = normalise_to_b0_and_temin(img,gradechoinv)
%normalise diffusion-relaxometry images so that the mean of the b=0 images at the lowest echo time is
%equal to 1.
%
%inputs 
%img - the diffusion weighted image to normalise
%b - corresponding b-values
%te - corresponding echo times
%
%output
%normimg - the normalised image

% Author
% Paddy Slator (p.slator@ucl.ac.uk)

%extract b-values and echo times
b = gradechoinv(:,4);
te = gradechoinv(:,5);


%check that the number of b values matches the number of images 
if size(img,4)~=length(b)
   disp('can''t normalise dw image: number of b-values doesn''t match number of volumes')
   normimg=[];
   return
end

%add smallest number to the image to prevent divide by zero errors
img = img + eps;
%make sure it's a double
img = double(img);

%find b0 and min te image (or images)
b0_minte_index=find(b==0 & te==min(te));


%if more than one b0+min te image, take the mean of them
if length(b0_minte_index)>1
    b0_minte_image=mean(img(:,:,:,b0_minte_index),4);
elseif length(b0_minte_index)==1
    b0_minte_image=img(:,:,:,b0_minte_index);
else
    disp('can"t normalise dw image: no b0 volumes')
    normimg=[];
    return
end

normimg=zeros(size(img));

for i=1:length(b)
    normimg(:,:,:,i)=img(:,:,:,i)./b0_minte_image;
end

   


end