function plot_inspect_seg(output)


for i=1:output.options.nclus    
    plot_multidim_spectrum(output.ILT_output{i}.F,...
        output.ILT_output{i}.grid,...
        output.params);

    title(['Multidimensional spectrum for cluster ' num2str(i)])            
end


figure;
for i=1:length(output.MLroi)
    imagesc(output.MLroi{i}(:,:,1))
    colorbar
    
    title('Voxelwise clustering')
end



