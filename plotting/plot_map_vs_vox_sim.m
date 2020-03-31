function plot_map_vs_vox_sim(inspect_map_fit,vox_fit,ground_truth)

%unpack some useful inspect_map output
%grid to plot the spectra on
grid = getkernelgrid(inspect_map_fit.options.ILT);
%output spectra
Fcomp = inspect_map_fit.iter{end}{end}.Fcomp;
%inspect output voxelwise spectral weights
inspect_map_vfimg = inspect_map_fit.iter{end}{end}.imgweights;
%number of spectral components
ncomp = inspect_map_fit.options.ncomp;


%unpack some useful voxelwise output
vox_vfimg = vox_fit.vfimg;



% %plot the output spectra
% for i=1:ncomp
%     figure;hold on;
%     
%     plot_multidim_spectrum(Fcomp,grid,kernel)
%     %contour(grid{2},grid{1},Fcomp{i})
%     
%     plot(spectparams.T2(i),spectparams.D(i),'rx')
%     title(['Component ' num2str(i)])
%     legend({'InSpect fit','Ground Truth'})
%     set(gca, 'YScale', 'log');
%     xlabel('T2* (s)')
%     ylabel('ADC (mm^2/s)')
% end

%plot all the maps in a big figure
figure; hold on;
subfigdim = [3 ncomp];
plotorder = 1:subfigdim(2);


for i=1:ncomp
    %plot the inspect_map maps
    %subtightplot(subfigdim(1),subfigdim(2),i)
    subplot(subfigdim(1),subfigdim(2),i)
    imagesc(inspect_map_vfimg{1}(:,:,1,plotorder(i)))
    
    colorbar
    %caxis([0 cmax(i)])
    if i==1
        ylabel('InSpect Maps')
    end

    %plot the ground truth maps
    %subtightplot(subfigdim(1),subfigdim(2),i+ncomp)
    subplot(subfigdim(1),subfigdim(2),i + ncomp)
    imagesc(ground_truth(:,:,1,i))
    colorbar;
    if i==1
        ylabel('Ground truth volume fraction')
    end

    %subtightplot(subfigdim(1),subfigdim(2),i + 2 * ncomp)
    subplot(subfigdim(1),subfigdim(2),i + 2 * ncomp)
    imagesc(vox_vfimg{1}(:,:,:,plotorder(i)))
    colorbar
    %caxis([0 cmax(i)])

    if i==1
        ylabel('Voxelwise Maps')
    end

end






%hard-coded!
set(gcf,'Position',[-56 1011 1261 677])


