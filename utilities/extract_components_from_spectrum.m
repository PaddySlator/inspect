function [Fcomp,Fcompvec] = extract_components_from_spectrum(F)
%identify the separated peaks in a spectrum

disp('want to extend this so it takes number of peaks as input - then adaptively change thresh until you get this many peaks')

%lower sensitivity prevents two close peaks merging
sensitivity = 0.01;
binF = imbinarize(F,'adaptive','Sensitivity',sensitivity);       

%pixels/components are connected if their edges touch (no diagonals)
connectivity = 4;
connections = bwconncomp(binF,connectivity);  



npeaks = connections.NumObjects;

for i=1:connections.NumObjects
    peaks{i} = connections.PixelIdxList{i};
end

Fvec = F(:);

for i=1:npeaks
    Fcompvec{i} = zeros(size(Fvec));
    Fcompvec{i}(peaks{i}) = Fvec(peaks{i});

    Fcomp{i} = reshape(Fcompvec{i},size(F));    
end