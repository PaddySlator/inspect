function S = simulate_spectrum_signal(F,gradechoinv,options)


%get the number of measurements 
Nmeas = size(gradechoinv,1);



%build the matrix of kernels - i.e. the dictionary - Nmeas x w1*w2*w3
K = build_kernel_matrix(gradechoinv,options);


Fvec = F(:);

%sythesise the signal
S = K*Fvec;

if isfield(options,'SNR')
    v = randn(size(S))./SNR;% Gaussian random noise
    S = (S+v);% Noisy signal 
end






end