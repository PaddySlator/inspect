function [TI_shift, acq_shift_vector] = calculate_T1_kernal_from_superblocks(gradechoinvsuper)

%The inversion times (TI) are not straight lines

%get the number of each superblock 
SBlabels = gradechoinvsuper(:,8);

%get the acquisition order for all superblocks
SBorder = gradechoinvsuper(:,9);


%get the inversion times
TI = gradechoinvsuper(:,6); 


%assume that all the superblocks have the same acquisition strategy
i=1;


%get the acquisition order for a superblock - 
n = SBorder( SBlabels == i );

%get the shift in TIs between acquisitions

%assuming that it is the same for all acquisitions
acq_shift = TI(SBlabels == i & SBorder == 2 ) - TI(SBlabels == i  & SBorder == 1);
acq_shift = acq_shift(1);
%calculate the shifted/corrected TI's
TI_shift = TI(SBlabels == i) - (SBorder(SBlabels == i) - 1)* acq_shift;


figure;hold on;
plot(TI_shift,'o')
plot(TI(SBlabels == i),'x')

%make a vector of all shifts (not sure if this is feasible to fit into
%MERA?
%what the TI's would be if there was no shift 
TI_no_shift = repmat(TI(SBlabels == i & SBorder == min(SBorder)) , [max(SBorder) 1]);

acq_shift_vector = TI(SBlabels == i) - TI_no_shift;


TI_shift =  TI(SBlabels == i) - acq_shift_vector;


figure;hold on;
plot(TI_shift,'o')
plot(TI(SBlabels == i),'x')




end
