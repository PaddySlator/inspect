function sigma = estimate_sd(img,gradechoinv)

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

if onb && ~onTE && ~onTI %diffusion
    sigma = std(img(:, b == 0)');
elseif onb && onTE && ~onTI %T2-diffusion
    sigma = std(img(:, b == 0 & TE == min(TE))');
elseif onb && ~onTE && onTI %T1-diffusion
    sigma = std(img(:, b == 0 & TI > TIratio * TR)');
elseif ~onb && onTE && onTI %T1-T2
    sigma = std(img(:, TE == min(TE) & TI > TIratio * TR)');
elseif onb && onTE && onTI %T1-T2-diffusion        
    sigma = std(img(:, b == 0 & TE == min(TE) & TI > TIratio * TR)');    
end

%make sure that output is a column vector 
if ~iscolumn(sigma)
    sigma = sigma';
end
