function S = simulate_spectrum_signal(F,gradechoinv,options,SNR)


%get the number of measurements 
Nmeas = size(gradechoinv,1);

%get the sampling times 
%b-value
b1 = gradechoinv(:,4);
%echo time
b2 = gradechoinv(:,5);

%get the size of the grid on which the spectrum is defined
Nk1 = size(F,1);
Nk2 = size(F,2);

Nk = size(F);

%build the matrix of kernels - i.e. the dictionary - Nmeas x w1*w2*w3
kernels = {'expdecay','expdecayrecip'};
K = build_kernel_matrix(kernels,gradechoinv,options);


reg=1;
if reg %regularisation    
    %augment kernal matrix 
    H = speye(prod(Nk),prod(Nk));
    %H(end,end) = 0;  
    
    Kalpha = [K; options.alpha*H];
                   
end

Fvec = F(:);

%sythesise the signal
S = K*Fvec;

v = randn(size(S))./SNR;% Gaussian random noise
S = (S+v);% Noisy signal 

%unaugment data vector
S = S(1:Nmeas);





end