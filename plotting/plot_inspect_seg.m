function plot_inspect_seg(output)


%get the subplot dimensions
nclus = output.options.nclus;

nrow = floor(sqrt(nclus));
ncol = nclus/nrow;

figure;


for i=1:nclus    
    subplot(nrow,ncol,i);
    
    plot_multidim_spectrum(output.ILT_output{i}.F,...
        output.ILT_output{i}.grid,...
        output.kernel);

    title(['Cluster ' num2str(i) ' Spectrum'])            
end



for i=1:length(output.MLroi)
    figure;
    nslices = size(output.MLroi{i},3);
    
    for j=1:nslices
        subplot(1,nslices,j);
        imagesc(output.MLroi{i}(:,:,j))
        
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);

        title(['Clusters slice ' num2str(j)])        
        
        cb=colorbar;
        cb.Location = 'southoutside';        
        %cb.Position = cb.Position + 1e-100;
    end    
         
   
    %cb.Position = cb.Position + 1e-100;%this stops the image from being resized
    cmap = colormap('jet');
    
    set(gcf,'Position',[1 458 1280 240]) %hard-coded for my screen!
       
end



