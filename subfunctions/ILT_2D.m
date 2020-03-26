function output = ILT_2D(sig, gradechoinv,options)

%solve S = KF for F
% S is the signal  
% K is the matrix of kernals
% F is the spectrum  

%make the matrices

%get the number of measurements 
Nmeas = size(gradechoinv,1);

%get the sampling times 
%b-value
b1 = gradechoinv(:,4);
%echo time
b2 = gradechoinv(:,5);




%define the grid on which to estimate the diffusion/relaxometry spectrum
Nk1 = options.Nk1;
Nk2 = options.Nk2;


w1 = logspace(log10(options.mink1),log10(options.maxk1),Nk1);
w2 = linspace(options.mink2,options.maxk2, Nk2);








%build the signal matrix - Nmeas x 1 
if iscolumn(sig)
    S = sig;
else
    S = sig';
end

%define the kernals


%build the matrix of kernals - Nmeas x w1*w2*w3
K = zeros(Nmeas,Nk1*Nk2);

for meas = 1:Nmeas
    l=1;
    for j=1:Nk2
        for i=1:Nk1
            K(meas,l) = exp(-b1(meas)*w1(i)) * ...
                exp(-b2(meas)/w2(j));
            l=l+1;
        end
    end
end

if options.reg %regularisation
    
    %augment kernal matrix 
    H = eye(Nk1*Nk2,Nk1*Nk2);
    %H(end,end) = 0;
  
    
    Kalpha = [K; options.alpha*H];
   
    
    %augment data vector    
    S = [S; zeros(size(H,1),1)];
    
        
end


tic;
if options.reg
    [Fvec,RESNORM] = lsqnonneg(Kalpha,S);
else
    [Fvec,RESNORM] = lsqnonneg(K,S);
end
toc;

%reorder F
F = zeros(Nk1,Nk2);
l=1;
for j=1:Nk2
    for i=1:Nk1
        F(i,j) = Fvec(l);
        l=l+1;
    end
end


%calculate the residuals 
if options.reg
    res = Kalpha*Fvec - S;
else
    res = K*Fvec - S;
end

%estimate the SNR 
SNR = sum(Fvec)./std(res);

output.SNR = SNR;
output.res = res;


output.S = S;
%vector format
output.Fvec = Fvec;
output.F = F;
output.RESNORM = RESNORM;
output.w1 = w1;
output.w2 = w2;
%unaugmented
output.K = K;
if options.reg
    %augmented
    output.Kalpha = Kalpha;
end








end