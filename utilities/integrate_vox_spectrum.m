function Vf = integrate_vox_spectrum(F,grid,sROI)

%get the number of dimensions from the supplied boundaries 
ndim = length(sROI);
%get number of spectra points
nspectelem = length(F(:));

%get a vector, same length as the number of spectra points,
%with index indicating corresponding integration region 
gridvec = combvec(grid{:})';

l=1;

%get the number of sROIs
nSROI = size(sROI{1},1);
    
        
SROI = zeros(nspectelem,1);

for j=1:nSROI
    indexinSROI = zeros(nspectelem,ndim);
    for i=1:ndim %loop over dimensions          
        %find elements in the sROI bounds for this dimension
        indexinSROI(:,i) = sROI{i}(j,1) < gridvec(:,i) & gridvec(:,i) <  sROI{i}(j,2);        
    end    
    %elements in the sROI bounds for all dimensions are in the sROI    
    SROI(sum(indexinSROI,2) == ndim) = j;
end
   
Fvec = F(:);
Vf = zeros(1,nSROI);
for j=1:nSROI
    Vf(j) = sum(Fvec(SROI == j));   
end

%normalise
Vf = Vf./sum(Vf);   
