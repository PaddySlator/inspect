function plot_multidim_spectrum(F,grid,kernel)



params = GetKernelParameterStrings(kernel);

ndim = length(grid);

figure;hold on;

if ndim == 1 %1D spectrum
    plot(grid,F)
    xlabel(params{1})
elseif ndim == 2 %2D spectrum
    contour(grid{2},grid{1},F);
    xlabel(params{2})
    ylabel(params{1})   
elseif ndim >= 3 %3D and above - need to plot projections
    
    %plotting a 2D projection
    plotdim = 2;
    
    %get the dimensions to plot
    dim2plot = combntns(1:ndim,plotdim);
    
    %number of projections 
    nproj = size(dim2plot,1);
            
        
    for i=1:size(dim2plot,1)
        
        %get the dimensions that we "squeeze" for the plot by removing the
        %plotted dimensions
        dim2squeeze = 1:ndim;
        dim2squeeze(dim2plot(i,:))=[];
                
        Fproj=F;
        for j=dim2squeeze                
            %squeeze the spectrum in this dimension
            Fproj = sum(Fproj, j);
        end
        Fproj = squeeze(Fproj);
          
        %plot this 2D projection       
        subplot(1,nproj,i);axis square;
        contour(grid{dim2plot(i,2)},grid{dim2plot(i,1)},Fproj)        
        
        xlabel(params{dim2plot(i,2)})
        ylabel(params{dim2plot(i,1)})
        
    end        


    
    %hard-coded for my screen!    
    set(gcf, 'Position', [-304        1508        1920         277])
    
end

