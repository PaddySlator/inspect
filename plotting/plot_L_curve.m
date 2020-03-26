function plot_L_curve(sig,gradechoinv,options,alpha)

%alpha is the range of regularisation parameters to plot


for i=1:length(alpha)
    options.alpha = alpha(i);
    output{i} = ILT_3D(sig, gradechoinv,options);        
    
    %store variables for plot
    logRESNORM(i) = log(output{i}.RESNORM);
    logSOLNORM(i) = norm(output{i}.F(:));
    
    %plot
    plot_3D_spectrum(output{i})
end


%plot log residual norm against log solution norm 
figure;plot(logRESNORM,logSOLNORM);

figure;plot(alpha,logSOLNORM)



end