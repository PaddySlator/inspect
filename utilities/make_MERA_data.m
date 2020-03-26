function MERA_data = make_MERA_data(sig,gradechoinv)
%make a MERA structure for diffusion-relaxometry data from single measurement 
%INPUTS: sig - the diffusion-relaxometry signal 
%	gradechoinv - matrix with form [gx gy gz b te], 
% 	where gx, gy, gz are the gradient directions, b are the b-values, te the echo times
%
%OUTPUT: MERA_data - structure suitable for MERA data fitting with fields:
%	D - Nbvals x Nechotimes data matrix - columns have the same echo time, rows have the same b-value 
%	t - column vector of b-values
%   t2 - column vector of echo times

if size(gradechoinv,2) <= 6
    disp('2D case')
    te = unique(gradechoinv(:,5));
    %need all te's to have the same set of b-values - a requirement for MERA fitting function
    for j=1:length(te)
        bcheck(:,j) = gradechoinv(gradechoinv(:,5) == te(j),4);
    end
    %check if same b-values for each te
    if nnz(diff(bcheck')) == 0
        b = bcheck(:,1);
    else
        error('need the same set of b-values for each te for MERA fitting')
    end
    
    
    % make the data matrix - columns with the same echo time
    for j=1:length(te)
        MERA_data.D(:,j) =  sig(gradechoinv(:,5) == te(j))';
    end
    %b-values - rescale
    MERA_data.t = b * 10^-3;
    %echo times
    MERA_data.t2 = te;
    
else
    disp('3D case')
    
    %do the echo times and b-values first
    te = unique(gradechoinv(:,5));   
    
     %need all te's to have the same set of b-values - a requirement for MERA fitting function
    for j=1:length(te)
        bcheck(:,j) = gradechoinv(gradechoinv(:,5) == te(j),4);
    end
    %check if same b-values for each te
    if nnz(diff(bcheck')) == 0
        b = bcheck(:,1);
    else
        error('need the same set of b-values for each te for MERA fitting')
    end
    
    %b-values - rescale
    MERA_data.t = b * 10^-3;
    %echo times
    MERA_data.t2 = te; 
    
    
    %calculate the shifted inversion times 
    [TI_shifted, TI_shift_vector] = calculate_T1_kernal_from_superblocks(gradechoinv);
    
       
    %inversion times
    MERA_data.t3 = TI_shifted;  
    
    MERA_data.t3_shift_vector = TI_shift_vector;
    
    
    
    
    % make the data matrix - columns with the same echo time, rows with the
    % same b-values, depth with the same echo times
    MERA_data.D = zeros();
    for j=1:length(te)
        MERA_data.D(:,j) =  sig(gradechoinv(:,5) == te(j))';
    end
    
end







end