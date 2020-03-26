function S = add_noise(S,SNR,noisetype)
%add noise to MRI signal 

if SNR~=Inf
    
    if strcmpi(noisetype,'gaussian')
        v = randn(size(S))./SNR;% Gaussian random noise
        S = (S+v);% Noisy signal 
    end
    
    if strcmpi(noisetype,'rician')
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
