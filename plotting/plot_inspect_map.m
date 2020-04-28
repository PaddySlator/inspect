function plot_inspect_map(output,options)



ncomp = output.options.ncomp;

nrows = ncomp;

nimg = length(output.imgweights);






for i=1:nimg  
    figure;
    
    slicetoplot = find(squeeze(sum(sum(output.imgweights{i}(:,:,:,1))) ~=0));%don't plot blank slices
    
    %nslices = size(output.imgweights{i},3);%get the number of slices for this image
    nslices = length(slicetoplot);
    
    ncol = nslices + 1;
    
    %plot the spectra in the first column  

    for j=1:ncomp
        subplot(nrows,ncol,1 + ncol*(j-1));
        
        plot_multidim_spectrum(output.iter{end}{end}.Fcomp{j},...
            output.ILT_mean.grid,...
            output.options.kernel);
        
        ylabel(['Component ' num2str(j)])
        
       
        %crop any rows/columns that are all zeros in every slice
        for k=1:nslices
           %allvf = sum(output.imgweights{1}(:,:,k,:),4); 
           nonzerorowsbyslice(k,:) = (~all(output.imgweights{1}(:,:,slicetoplot(k)) == 0,2));
           nonzerocolsbyslice(k,:) = (~all(output.imgweights{1}(:,:,slicetoplot(k)) == 0,1));
        end
	       
        nonzerorows = sum(nonzerorowsbyslice,1) ~= 0 ;
        nonzerocols = sum(nonzerocolsbyslice,1) ~= 0 ;  
                
        
        for k=1:nslices            
            subplot(nrows,ncol,ncol*(j-1) + k + 1);                                
            
            window=5;
            imagesc(output.imgweights{i}(nonzerorows,nonzerocols,slicetoplot(k),j));
                        
            set(gca,'YTickLabel',[]);
            set(gca,'XTickLabel',[]);
                        
        end
        %cb=colorbar;
        %cb.Location = 'eastoutside';
        %cb.Position = cb.Position + [0 0 0 0]; %hard-coded!
    end
    
    caxis([0 1])
    cmap = colormap('parula');
    
    set(gcf,'Position',[60 28 1124 677]) %hard-coded for my screen!             
        
    if options.save
        print_to_formats([options.save_path options.dirname '/inspect_map_' num2str(options.ncomp) '_comp_' options.scan_names{i}],{'-depsc','fig','-dpng'})
    end
    
end
