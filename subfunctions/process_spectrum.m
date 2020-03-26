function [Fcomp,weights] = process_spectrum(F)
%takes a spectrum as input and identifies separated peaks 
%(i.e. "spectral components") and their respective volume fractions

%binarize the image
binF = imbinarize(F);        
%identify connected components in binary image
connections = bwconncomp(binF);    
%extract the position of each peak
for i=1:connections.NumObjects
    peakpos{i} = connections.PixelIdxList{i};
    %convert to coordinates
    %[x,y] = ind2sub(peakpos{i},size(F));
end

%need to make spectral compenents in vector form first
Fvec = F(:);
Fcompvec = cell(1,connections.NumObjects);    
Fcomp = cell(1,connections.NumObjects);    

weights = zeros(1,connections.NumObjects);    

  
for i=1:connections.NumObjects
   %peak pos are in ind form
   Fcompvec{i} = zeros(size(F(:))); 
   
   Fcompvec{i}(peakpos{i}) = Fvec(peakpos{i});
   
   Fcomp{i} = reshape(Fcompvec{i}, size(F));
   
   %get the weight of this spectral component 
   weights(i) = sum(Fcomp{i}(:))/sum(F(:));
end
   
%normalise the weights
weights = weights./sum(weights);
    
    


end



