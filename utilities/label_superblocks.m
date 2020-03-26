function gradechoinvsuper = label_superblocks(gradechoinv)
%label ZEBRA "superblocks"


gradechoinvsuper = gradechoinv;


superblockstarts = find(gradechoinv(:,6) == min(gradechoinv(:,6)));

%get the number of slices by finding the gap between lowest TI's
Ns = diff(superblockstarts(1:2));

%get the number of superblocks
Nb = length(superblockstarts);

%label each superblock
sliceorder = repmat(1:Nb,[Ns 1]);
sliceorder = sliceorder(:);

gradechoinvsuper(:,8) = sliceorder;


%label the acquisition order (I'm counting each inversion pulse as a new "acquisition") for each superblock 

%find the number of slices for each inversion pulse
inversionpulses = find( gradechoinv(1:Ns-1,6) > gradechoinv(2:Ns,6) );
%number of slices read out for each inversion pulse
Nsinv = inversionpulses(1);
%number of inversion pulses for each superblock
Ninv = Ns / Nsinv;

acquisitionorder = repmat(1:Ninv,[Nsinv 1]);
acquisitionorder = acquisitionorder(:);

%label the acquisition order for all superblocks
gradechoinvsuper(:,9) =  repmat(acquisitionorder,[Nb 1]);

end











