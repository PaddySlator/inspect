


%three peak model
D = [0.01 0.01 0.001];
T2 = [200 50 200];

D = [1/200 1/100 1/100];
T2 = [200 100 200];

%volume fractions
f1 = 0.3;f2 = 0.3;
F = [f1 f2 1-f1-f2];


%sampling times - e.g. b-values and echo times 
b=50:1:250;
te = 50:1:250;


% %one peak model
% D=[0.05 0.05 0.05];
% T2=[0.35 0.35 0.35];
% 
% b = .001*(1:128)';
% te = logspace(log10(0.02),log10(5),64)'; 


%simulate the experiment
S = zeros(length(b), length(te));

for i=1:length(b)
    for j=1:length(te)
        S(i,j) = F(1)*( exp( -te(j)/T2(1) ) * exp( -b(i)*D(1) ) ) ...
            + F(2)*( exp( -te(j)/T2(2) ) * exp( -b(i)*D(2) ) ) ...
            + F(3)*( exp( -te(j)/T2(3) ) * exp( -b(i)*D(3) ) );
    end
end


SNR = 50000;

v = randn(size(S))./SNR;% Gaussian random noise
S = (S+v);% Noisy signal 
data.D = S;
data.t = b;
data.t2 = te;


%% Calling MERA
% profile -memory on
clear fitting
close all
fitting.twoD = 'y';
fitting.threeD = 'n';
analysis.interactive = 'n';
fitting.regtyp = 'me';
fitting.rangeT = [50, 250];
fitting.rangeT2 = [50, 250];
fitting.numbergauss = 2;
fitting.regadj = 'manual';
analysis.graph = 'y';
fitting.regweight = 0.002;
fitting.numberT = 40;
fitting.numberT2 = 40;
analysis.extract = 'auto';

[MERA_out2D,MERA_fitting_out]=MERA(data,fitting,analysis);



    
    