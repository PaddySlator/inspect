function plot_inspect_seg(output)


for i=1:output.options.nclus    
    plot_multidim_spectrum(output.ILT_output{i}.F,...
        output.ILT_output{i}.grid,...
        output.kernel);

    title(['Multidimensional spectrum for cluster ' num2str(i)])            
end


figure;
for i=1:length(output.MLroi)
    
    nslices = size(output.MLroi{i},3);
    
    for j=1:nslices
        imagesc(output.MLroi{i}(:,:,1))
        colorbar
    end
    
    title('Voxelwise clustering')
end



