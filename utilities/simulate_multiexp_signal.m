function S = simulate_multiexp_signal(params,gradechoinv,options)

%params, kernal are cells with 


%simulate an experiment using the gradechoinv table


%get the sampling times 
b = gradechoinv(:,4);
te = gradechoinv(:,5);

%get the SNR and noise type
SNR = options.SNR;
noisetype = options.noisetype;

%define the diffusion/relaxometry constants 
%three compartment model
D = params.D; 
T2 = params.T2;
f = params.f;


S = zeros(length(gradechoinv),1);


for i=1:length(S)
    S(i)=0;
    for j = 1:length(f)
        S(i) = S(i) + ...
            f(j) * (...
                exp(-b(i)*D(j)) * ...
                exp(-te(i)/T2(j))...
                );
                %abs(1 - 2*exp(-(b3(i) + b2(i))/T1(j)) + exp(-TR/b3(i)))...
                %);
    end
end


if SNR~=Inf
    
    if strcmp(noisetype,'gaussian')
        v = randn(size(S))./SNR;% Gaussian random noise
        S = (S+v);% Noisy signal 
    end
    
    if strcmp(noisetype,'rician')
        %Rician noise
        sigma = 1/SNR;
        realnoise = randn(size(S))*sigma;
        imagnoise = randn(size(S))*sigma;
        
        %noisy complex data
        Scomplex = S + realnoise + 1i * imagnoise;
        
        %rician distributed magnitude data 
        S = abs(Scomplex);
        
        
    end
end









end