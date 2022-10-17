function S = KernelDD(params,gradechoinv)
%kernel for double diffusion encoding experiment
%This is a bit of a hack where the second b-value is taken from the 8th column 
%of the gradechoinv file - or maybe it's ok to keep it like this and define
%the 8th column as b2?

%equation is 
%this is for a fixed mixing time

d1 = params(1);
d2 = params(2);

%put the b2-values into a separate gradechoinv
gradechoinv2 = gradechoinv;
gradechoinv2(:,4) = gradechoinv(:,8); 

S = KernelD(d1,gradechoinv) .* KernelD(d2,gradechoinv2);

