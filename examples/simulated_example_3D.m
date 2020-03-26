


%three peak model
D = [0.001 0.001 0.003];
D=1./D;
D = [300 5000 5000];
T2 = [100 200 100];
T1 = [1000 2500 2500];


D = [100 200 200];
T2 = [100 200 100];
T1 = [200 100 100];

%D = [175 175 175];
%T2 = [175 175 175];
%T1 = [175 175 175];


%volume fractions
f1 = 0.3;f2 = 0.3;f3 = 0.2;
%F = [f1 f2 f3 1-f1-f2-f3];
F= [f1 f2 1-f1-f2];

%sampling times - e.g. b-values and echo times 
b=50:10:250;
te = 50:10:250;
ti = 50:10:250;

% zebra acquisition parameters
%b=[0 100 1000 2500];
%te = [57 81 171 228 285];
%ti = 500:500:5000;

% %one peak model
% D=[0.05 0.05 0.05];
% T2=[0.35 0.35 0.35];
% 
% b = .001*(1:128)';
% te = logspace(log10(0.02),log10(5),64)'; 


%simulate the experiment
S = zeros(length(b), length(te),length(ti));

for i=1:length(b)
    for j=1:length(te)
        for k=1:length(ti)
%             S(i,j,k) = F(1)*( exp( -te(j)/T2(1) ) * exp( -b(i)*D(1) ) * (1 - 2*exp(-ti(k)/T1(1)))) ...
%                 + F(2)*( exp( -te(j)/T2(2) ) * exp( -b(i)*D(2) ) * (1 - 2*exp(-ti(k)/T1(2)))) ...
%                 + F(3)*( exp( -te(j)/T2(3) ) * exp( -b(i)*D(3) ) * (1 - 2*exp(-ti(k)/T1(3)))) ...
%                 + F(4)*( exp( -te(j)/T2(4) ) * exp( -b(i)*D(4) ) * (1 - 2*exp(-ti(k)/T1(4))));

           % S(i,j,k) = F(1)*( exp( -te(j)/T2(1) ) * exp( -b(i)*D(1) ) * exp(-ti(k)/T1(1))) ...
           %      + F(2)*( exp( -te(j)/T2(2) ) * exp( -b(i)*D(2) ) * exp(-ti(k)/T1(2))) ...
           %      + F(3)*( exp( -te(j)/T2(3) ) * exp( -b(i)*D(3) ) * exp(-ti(k)/T1(3))); ...
           %+ F(4)*( exp( -te(j)/T2(4) ) * exp( -b(i)*D(4) ) * exp(-ti(k)/T1(4)));
           
           
           S(i,j,k) = F(1)*( exp( -te(j)/T2(1) ) * exp( -b(i)/D(1) ) * exp(-ti(k)/T1(1))) ...
               + F(2)*( exp( -te(j)/T2(2) ) * exp( -b(i)/D(2) ) * exp(-ti(k)/T1(2))) ...
               + F(3)*( exp( -te(j)/T2(3) ) * exp( -b(i)/D(3) ) * exp(-ti(k)/T1(3))); ...
        end
    end
end


SNR = 50000;

v = randn(size(S))./SNR;% Gaussian random noise
S = (S+v);% Noisy signal 
data.D = S;
data.t = b;
data.t2 = te;
data.t3 = ti;



% %% Calling MERA
% % profile -memory on
% clear fitting
% close all
% fitting.twoD = 'n';
% fitting.threeD = 'y';
% analysis.interactive = 'n';
% fitting.regtyp = 'me';
% fitting.rangeT = [50, 250];
% fitting.rangeT2 = [50,250];
% fitting.rangeT3 = [50,250];
% fitting.numbergauss = 3;
% fitting.regadj = 'manual';
% analysis.graph = 'y';
% fitting.regweight = 0.002;
% fitting.numberT = 10;
% fitting.numberT2 = 10;
% fitting.numberT3 = 10;
% analysis.extract = 'auto';
% fitting.nnlscode = 'lsqnonneg';
% fitting.nnlscode = 'nnlsmex';
% %fitting.nnlscode = 'fnnls';
% 
% fitting.theta_vector = 180;
% 
% 
% [output,MERA_fitting_out]=MERA_3D(data,fitting,analysis);
% 
% 
% 
% 
% proj = squeeze(sum(output.S,1));
% figure;hold on;
% contour(output.T3,output.T2,proj)
% xlabel('T3')
% ylabel('T2')
% %add the simulation peaks
% for i=1:length(T1)
%     plot(T1(i),T2(i),'rx')
% end
% 
% 
% proj = squeeze(sum(output.S,2));
% figure;hold on;
% contour(output.T3,output.T,proj)
% xlabel('T3')
% ylabel('T')
% for i=1:length(T1)
%     plot(T1(i),D(i),'rx')
% end
% 
% 
% 
% proj = squeeze(sum(output.S,3));
% figure;hold on;
% contour(output.T2,output.T,proj)
% xlabel('T2')
% ylabel('T')
% for i=1:length(T1)
%     plot(T2(i),D(i),'rx')
% end
% 
% 
% 
% %plot the synthetic data, S
% % proj = squeeze(sum(S,1));
% % figure;
% % surf(ti,te,proj)
% % xlabel('ti')
% % ylabel('te')
% % 
% % proj = squeeze(sum(S,2));
% % figure;
% % surf(ti,b,proj)
% % xlabel('ti')
% % ylabel('b')
% % 
% % 
% % proj = squeeze(sum(S,3));
% % figure;
% % surf(te,b,proj)
% % xlabel('te')
% % ylabel('b')




%% try with my new code!
datadir = '/Users/paddyslator/Dropbox/placentaJhu/Data_other/t1diff/zebra101';

gradechoinv = load([datadir '/gradechoinvTR.txt']);

%simulate an experiment using the gradechoinv table
%get the sampling times 
b1 = gradechoinv(:,4);
b2 = gradechoinv(:,5);
b3 = gradechoinv(:,6);

%define the diffusion/relaxometry constants 
%three compartment model
D = [0.001 0.0003 0.0003] ;
T2 = [100 150 150];
T1 = [1000 3000 1000];
%volume fractions
f1 = 0.3;f2 = 0.3;f3 = 0.2;
%F = [f1 f2 f3 1-f1-f2-f3];
f = [f1 f2 1-f1-f2];

TR = 10000;

S = zeros(length(gradechoinv),1);


for i=1:length(S)
    S(i)=0;
    for j = 1:length(f)
        S(i) = S(i) + ...
            f(j) * (...
                exp(-b1(i)*D(j)) * ...
                exp(-b2(i)/T2(j)) * ...
                abs(1 - 2*exp(-(b3(i) + b2(i))/T1(j)) + exp(-TR/b3(i)))...
                );
    end
end

SNR = 50000;

v = randn(size(S))./SNR;% Gaussian random noise
S = (S+v);% Noisy signal 


ILT_options.Nk1=30;
ILT_options.Nk2=30;
ILT_options.Nk3=30;

ILT_options.mink1 = 0.00001;
ILT_options.maxk1 = 0.0015;
ILT_options.mink2 = 50;
ILT_options.maxk2 = 200;
ILT_options.mink3 = 500;
ILT_options.maxk3 = 3500;


ILT_options.alpha = 0.01;


output=ILT_3D(S,gradechoinv,ILT_options);



%plot projections
proj = squeeze(sum(output.F,1));
figure;hold on;
contour(output.w3,output.w2,proj)
xlabel('T1')
ylabel('T2')
%simulation peaks
for i=1:length(T1)
    plot(T1(i),T2(i),'rx')
end


proj = squeeze(sum(output.F,2));
figure;hold on;
contour(output.w3,output.w1,proj)
xlabel('T1')
ylabel('D')
%simulation peaks
for i=1:length(T1)
    plot(T1(i),D(i),'rx')
end

proj = squeeze(sum(output.F,3));
figure;hold on;
contour(output.w2,output.w1,proj)
xlabel('T2')
ylabel('D')
%simulation peaks
for i=1:length(T1)
    plot(T2(i),D(i),'rx')
end




