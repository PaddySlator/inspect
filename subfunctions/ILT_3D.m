function output = ILT_3D(sig, gradechoinv,options)

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
%inversion time
b3 = gradechoinv(:,6);
%repetition time
TR = gradechoinv(1,7);

%define the grid on which to estimate the diffusion/relaxometry spectrum
Nk1 = options.Nk1;
Nk2 = options.Nk2;
Nk3 = options.Nk3;

w1 = logspace(log10(options.mink1),log10(options.maxk1),Nk1);
w2 = linspace(options.mink2,options.maxk2, Nk2);
w3 = logspace(log10(options.mink3),log10(options.maxk3),Nk3);






%build the signal matrix - Nmeas x 1 
if iscolumn(sig)
    S = sig;
else
    S = sig';
end

%define the kernals


%build the matrix of kernals - Nmeas x w1*w2*w3
K = zeros(Nmeas,Nk1*Nk2*Nk3);


for meas = 1:Nmeas
    l=1;
    for k=1:Nk3
        for j=1:Nk2
            for i=1:Nk1
                K(meas,l) = exp(-b1(meas)*w1(i)) * ...
                    exp(-b2(meas)/w2(j)) * ...
                    abs(1 - 2*exp(-(b3(meas) + b2(meas))/w3(k)) + exp(-TR/w3(k)) );
                    ...abs(1 - 2*exp(-(b3(meas))/w3(k)) + exp(-TR/w3(k)) );
                    %abs(1 - 2*exp(-b3(meas)/w3(k)) );
                l=l+1;
            end
        end
    end
end

reg=1;
if reg %regularisation
    
    %augment kernal matrix 
    H = speye(Nk1*Nk2*Nk3,Nk1*Nk2*Nk3);
    %H(end,end) = 0;
  
    
    K = [K; options.alpha*H];
   
    
    %augment data vector    
    S = [S; zeros(size(H,1),1)];
    
        
end

tic;
[F,RESNORM] = lsqnonneg(K,S);
toc;

%reorder F
l=1;
for k=1:Nk3
    for j=1:Nk2
        for i=1:Nk1
            Fmat(i,j,k) = F(l);
            l=l+1;
        end
    end
end
F=Fmat;


output.S = S;
output.F = F;
output.RESNORM = RESNORM;
output.w1 = w1;
output.w2 = w2;
output.w3 = w3;









end