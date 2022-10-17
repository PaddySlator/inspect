function [sigma, S0] = estimate_sd(img,gradechoinv)

%estimate the standard deviation of image noise


%get acquisition parameters from the gradechoinv file
[~, b, TE, TI, TR] = unpack_gradechoinv(gradechoinv);


%check which acquisition parameters are included in this image - if the
%parameter isn't included the gradechoinv is all NaN's
onb = ~any(isnan(b));
onTE = ~any(isnan(TE));
onTI = ~any(isnan(TI));

%estimate the standard deviation by taking the standard deviation of the
%volumes with highest signal. We use  b=0 for diffusion, lowest TE for
%T2-decay, all volumes with TI > TIratio * TR for T1 inversion recovery.

%set this - visual inspection suggests 0.5 is good for ZEBRA data
TIratio = 0.5;

try
    if onb && ~onTE && ~onTI %diffusion        
        S0 = mean(img(:, b == min(b))');
        sigma = std(img(:, b == min(b))');        
    elseif onb && onTE && ~onTI %T2-diffusion
        S0 = mean(img(:, b == min(b) & TE == min(TE))');
        sigma = std(img(:, b == min(b) & TE == min(TE))');
    elseif onb && ~onTE && onTI %T1-diffusion
        S0 = mean(img(:, b == min(b) & TI > TIratio * TR)');
        sigma = std(img(:, b == min(b) & TI > TIratio * TR)');
    elseif ~onb && onTE && onTI %T1-T2
        S0 = mean(img(:, TE == min(TE) & TI > TIratio * TR)');
        sigma = std(img(:, TE == min(TE) & TI > TIratio * TR)');
    elseif onb && onTE && onTI %T1-T2-diffusion
        S0 = mean(img(:, b == min(b) & TE == min(TE) & TI > TIratio * TR)');
        sigma = std(img(:, b == min(b) & TE == min(TE) & TI > TIratio * TR)');
    end
catch % if there are not enough to take the standard deviation
    sigma = NaN;
end

%make sure that output is a column vector 
if ~iscolumn(sigma)
    sigma = sigma';
end
if ~iscolumn(S0)
    S0 = S0';
end