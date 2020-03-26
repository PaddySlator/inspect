datadir = '/Users/paddyslator/Dropbox/placentaJhu/Data_other/t1diff/';
datadir = '/Users/paddyslator/Dropbox/challenge/';


scans = {'zebra101','zebra102','zebra103'};
scans = {'zebra103'};

scans = {'cdmri0011','cdmri0012','cdmri0013','cdmri0014','cdmri0015','cdmri0016'};
scans = {'cdmri0011','cdmri0012'};
    
saveon = 0;

savedir = '~/Desktop/';

fig_formats = {'-dpdf','-dpng','fig'};

loadon = 1;

if loadon
    %load the scans
    clear input gradechoinv mask
    
    for i=1:length(scans)
        %load image
        imgfile = dir([datadir scans{i} '/MB_Re_t.nii*']);
        input.(scans{i}) = load_nii([datadir scans{i} '/' imgfile.name]);
        input.(scans{i}).img = double(input.(scans{i}).img) + eps;
        %gradechoinvfile
        %gradechoinvfile = dir([datadir scans{i} '/*paddy_correction*.csv']);
        %gradechoinvfile = dir('~/Desktop/*paddy_correction*.csv');
        gradechoinvfile = dir([datadir '/*parameters.txt*']);
        
        %gradechoinv.(scans{i}) = csvread([datadir scans{i} '/' gradechoinvfile.name]);
        gradechoinv.(scans{i}) = importdata([datadir '/' gradechoinvfile.name],' ');
        
        %these are actually in gradinvecho order so swap!
        gradechoinvtemp = gradechoinv.(scans{i});
        gradechoinv.(scans{i})(:,5) = gradechoinvtemp(:,6);
        gradechoinv.(scans{i})(:,6) = gradechoinvtemp(:,5);
        
        %mask
        mask.(scans{i}) = load_nii([datadir scans{i} '/MB_Re_t_mask.nii.gz']);
        mask.(scans{i}) = logical(mask.(scans{i}).img);
    end
    
    %add the TRs as a final column
    TR = 7500;
    for i=1:length(scans)
        gradechoinv.(scans{i})(:,7) = TR;
    end
    
end

%% normalise the images
for i=1:length(scans)
    norminput.(scans{i}) = input.(scans{i});
    norminput.(scans{i}).img = normalise_to_b0_temin_timax(input.(scans{i}).img,gradechoinv.(scans{i}));
    
end



%% plot the gradechoinv
figure; hold on
i=1;

%plot first ones for legend
for j=[21 5 1]
    if gradechoinv.(scans{i})(j,4) == 0
        plot(gradechoinv.(scans{i})(j,5),gradechoinv.(scans{i})(j,6),'bo')
    elseif gradechoinv.(scans{i})(j,4) == 2000
        plot(gradechoinv.(scans{i})(j,5),gradechoinv.(scans{i})(j,6),'rx')
    elseif gradechoinv.(scans{i})(j,4) == 3000
        plot(gradechoinv.(scans{i})(j,5),gradechoinv.(scans{i})(j,6),'g*')
    end
end


for j=1:size(gradechoinv.(scans{i}),1)
    if gradechoinv.(scans{i})(j,4) == 0
        plot(gradechoinv.(scans{i})(j,5),gradechoinv.(scans{i})(j,6),'bo')
    elseif gradechoinv.(scans{i})(j,4) == 2000
        plot(gradechoinv.(scans{i})(j,5),gradechoinv.(scans{i})(j,6),'rx')
    elseif gradechoinv.(scans{i})(j,4) == 3000
        plot(gradechoinv.(scans{i})(j,5),gradechoinv.(scans{i})(j,6),'g*')
    end
end

ylabel('Inversion time (TI)','FontSize',16)
xlabel('Echo time (TI)','FontSize',16)
legend({'b=0 s/mm^2','b=2000 s/mm^2','b=3000 s/mm^2'},'Location','NorthOutside','Orientation','horizontal','FontSize',16)

if saveon
    print_to_formats('/Users/paddyslator/Documents/OtherDocs/abstracts/gradechoinv',fig_formats)   
end





%% calculate the mean signal in the whole brain mask and in each slice
for i=1:length(scans)
    mean_brain_sig.(scans{i}) = calculate_mean_signal(input.(scans{i}).img,mask.(scans{i}));
    
    for j=1:input.(scans{i}).hdr.dime.dim(4)
        slicemask = zeros(size(mask.(scans{i})));
        slicemask(:,:,j) = mask.(scans{i})(:,:,j);
        mean_slice_sig.(scans{i}){j} = calculate_mean_signal(input.(scans{i}).img,slicemask);
    end
end

%% ilt options

ILT_options.Nk1=20;
ILT_options.Nk2=20;
ILT_options.Nk3=20;

ILT_options.Nk = [20 20 20];

%using camino limits now!
ILT_options.mink1 = 2 * 10^-4;
ILT_options.maxk1 = 0.01;
ILT_options.mink2 = 5;
ILT_options.maxk2 = 200;
ILT_options.mink3 = 250;
ILT_options.maxk3 = 4000;

ILT_options.mink = [0.00005*10^-6 5 250];
ILT_options.maxk = [0.01*10^-6 200 4000];


ILT_options.loggrid = [1 0 1];


ILT_options.reg = 1;
ILT_options.alpha = 0.01;

ILT_options.kernel = 'BallT2T1';


%% do ilt on each slice
for i=1:length(scans)
    for j=1:input.(scans{i}).hdr.dime.dim(4)    
        %slice_spect.(scans{i}){j}=ILT_3D(mean_slice_sig.(scans{i}){j},gradechoinv.(scans{i}),ILT_options);        
        slice_spect.(scans{i}){j}=ILT(mean_slice_sig.(scans{i}){j},gradechoinv.(scans{i}),ILT_options);
    end
end

%% plot some slices
zslice = [12 32 31];
i=1;j=12;
%figure;imagesc(input.(scans{i}).img(:,:,j,6))
plot_3D_spectrum(slice_spect.(scans{i}){j})
   
i=2;j=32;
%figure;imagesc(input.(scans{i}).img(:,:,j,6))
plot_3D_spectrum(slice_spect.(scans{i}){j})

i=3;j=31;
%figure;imagesc(input.(scans{i}).img(:,:,j,6))
plot_3D_spectrum(slice_spect.(scans{i}){j})

%% do ilt on whole brain masks
for i=1:length(scans)
    brain_spect.(scans{i})=ILT(mean_brain_sig.(scans{i}),gradechoinv.(scans{i}),ILT_options);
end

%% plot them
for i=1:length(scans)
   plot_3D_spectrum(brain_spect.(scans{i}))
end




%%
plot_3D_spectrum(slice_spect.(scans{i}){j})






%% run inspect segmentation version
inspect_options.ILT = ILT_options;

inspect_options.EM.n_clusters = 4;
inspect_options.EM.n_steps = 5;
inspect_options.EM.init = 'random';

inspect_options.save=0;
%inspect_options.save_path = paperpath;

for i=1:length(scans)
    [inspectoutput.(scans{i}), inspectoutputsummary.(scans{i})] = inspect_seg(input.(scans{i}).img,...
                    gradechoinv.(scans{i}),...
                    mask.(scans{i}),...
                    inspect_options);  
end



return 




% %% Phase correction
% 
% f = exp(-(-10:10).^2 / (2*4.^2));
% f = f./sum(f);
% ph = zeros(Nx,Ny,Nz,Nd);
% for d=1:Nd
%     for z=1:Nz
%         ph(:,:,z,d) = angle(conv2(f, f, resortedTI(:,:,z,d), 'same'));
% %                     if (mean(cos(p(mask(:,:,z)))) >= 0)
% %                         ph(:,:,z,d,t) = p;
% %                     else
% %                         ph(:,:,z,d,t) = p + pi;
% %                     end
%     end
% end
% 
% resortedTI = resortedTI .* exp(-1i.*ph);
% 
% inputcor = input;
% inputcor.img = resortedTI;


% Reject shortest TI




%save([datadir remove_ext_from_nifti(imgfile) '_' remove_ext_from_nifti(maskfile) '_T1-T2star-ADC_spectrum'],'output')


% plot projections
figure;hold on;

i=1;
subtightplot(1,3,i,0.05,0.1, 0.05)
proj = squeeze(sum(output.F,1));
contour(output.w3,output.w2,proj)
xlabel('T1')
ylabel('T2*')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

i=2;
subtightplot(1,3,i,0.05,0.1, 0.05)
proj = squeeze(sum(output.F,2));
contour(output.w3,output.w1,proj)
xlabel('T1')
ylabel('D')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

i=3;
subtightplot(1,3,i,0.05,0.1, 0.05)
proj = squeeze(sum(output.F,3));
contour(output.w2,output.w1,proj)
xlabel('T2*')
ylabel('D')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

% %plot scatter
% scatterpoints = spectrum_to_scatter(output);
% figure;
% plot3(scatterpoints(:,3),scatterpoints(:,2),scatterpoints(:,1),'o')
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% set(gca, 'ZScale', 'log')
% 
% xlabel('T1')
% ylabel('T2*')
% zlabel('D')






%% fit ILT voxelwise


i=1;
zslice=32;
nvox = nnz(mask.(scans{i})(:,:,zslice));
ILT_output=cell(nvox,1);

ILT_options.Nk1=20;
ILT_options.Nk2=20;
ILT_options.Nk3=20;

ILT_options.mink1 = 0.00005;
ILT_options.maxk1 = 0.01;
ILT_options.mink2 = 5;
ILT_options.maxk2 = 200;
ILT_options.mink3 = 250;
ILT_options.maxk3 = 4000;

ILT_options.alpha = 0.01;

clear ILT_output

dims = input.(scans{i}).hdr.dime.dim(2:5);

Nx=dims(1); Ny=dims(2); Nz=dims(3);
Nd=dims(4); 

l=1;
start = tic;
for z=zslice
    for y=1:Ny
        for x=1:Nx        
            if mask.(scans{i})(x,y,z)
                data = abs(squeeze(input.(scans{i}).img(x,y,z,:,:)));                              
                
                ILT_output{l}=ILT_3D(data,gradechoinv.(scans{i}),ILT_options);
                ILT_output{l}.voxel = [x y z];
                disp(['voxel ' num2str(l) ' of ' num2str(nvox) ]);
                l=l+1;                               
            end
        end
    end
end
runtime = toc(start);


%save([datadir 'voxelwise_spect_5_nov' remove_ext_from_nifti(maskfile) ],'ILT_output')
if saveon
%    save([savedir scans{i} '/voxelwise_spectrum_' num2str(zslice)],'ILT_output')
     save([savedir '/voxelwise_spectrum_' num2str(zslice)],'ILT_output')

end

%% fit inspect segmentation version

%put images and masks into cells 
for i=1:length(scans)
    img{i} = input.(scans{i}).img;
    masks{i} = mask.(scans{i});
end

%they all have the same gradechoinv
gradechoinvall = gradechoinv.(scans{1});
        
        
inspect_options.ILT = ILT_options;
inspect_options.EM.n_clusters = 10;
inspect_options.EM.n_steps = 5;
inspect_options.EM.init = 'random';
inspect_options.EM.init = 'kmeans';

inspect_options.save = 0;

inspectsegzebra = inspect_seg(img,gradechoinvall,masks,inspect_options);





%%



%w1
ADC_bounds = [ILT_options.mink1 2*10^-4;
    ILT_options.mink1 ILT_options.maxk1;
    ILT_options.mink1 ILT_options.maxk1;
    ILT_options.mink1 ILT_options.maxk1
    10^-3 ILT_options.maxk1];

%w2
T2_bounds = [ILT_options.mink2 50;
    ILT_options.mink2 ILT_options.maxk2;
    ILT_options.mink2 ILT_options.maxk2; 
    ILT_options.mink2 ILT_options.maxk2;
    80 ILT_options.maxk2];

%w3
T1_bounds = [ILT_options.mink3 750; 
    750 1200; 
    1200 1850;
    1850 2500
    2500 5000];

zslice = [12 32 31];

%three compartments - white matter, grey matter, CSF
%w1
ADC_bounds = [ILT_options.mink1 10^-3;
    10^-3 ILT_options.maxk1;
    ILT_options.mink1 ILT_options.maxk1;
    ILT_options.mink1 ILT_options.maxk1;
    10^-3 ILT_options.maxk1];

%w2
T2_bounds = [ILT_options.mink2 ILT_options.maxk2
    ILT_options.mink2 ILT_options.maxk2;
    ILT_options.mink2 ILT_options.maxk2;
    ILT_options.mink2 ILT_options.maxk2; 
    70 ILT_options.maxk2];

%w3
T1_bounds = [ILT_options.mink3 650; 
    ILT_options.mink3 650; 
    650 1150;
    1150 2000; 
    2000 ILT_options.maxk3];

saveon = 0;

compartment_titles = {'White matter','White matter', 'Grey matter','CSF'};

for i=1:length(scans)
    %plot boundaries
    integrate_spectrum_3D(slice_spect.(scans{i}){zslice(i)},ADC_bounds,T2_bounds,T1_bounds,1)
        
    set(gcf,'Position',[ 82         328        1075         377])
    
    if saveon
        print_to_formats(['/Users/paddyslator/Documents/OtherDocs/abstracts/brain_spectrum_' scans{i} ],fig_formats)
    end
    
end


%do a separate legend for this
figure;hold on;
colors= hsv(length(compartment_titles));
for i=1:length(compartment_titles)
   plot(0,0,'color',colors(i,:),'Linewidth',3) 
end
legend(compartment_titles,'Fontsize',20,'NumColumns',3)


if saveon
   print_to_formats('/Users/paddyslator/Documents/OtherDocs/abstracts/brain_spectrum_legend',fig_formats)   
end





%% integrate over these boundaries
%number of compartments to integrate over
Ncomp = size(ADC_bounds,1);

saveon = 0;


%loadVfmap zebra101

zslice = 32;

clear voxel_spect Vfmap

for i=1:length(scans)
    %load([savedir scans{i} '/voxelwise_spectrum_' num2str(zslice(i)) '.mat'])
    load([savedir '/voxelwise_spectrum_' num2str(zslice(i)) '.mat'])
    voxel_spect{i} = ILT_output;
    
    Vfmap{i} = zeros([dims(1:3) size(ADC_bounds,1)]);
    
    %
    nvox=length(ILT_output);
    
    for j=1:nvox
        Vf = integrate_spectrum_3D(voxel_spect{i}{j},ADC_bounds,T2_bounds,T1_bounds,0);
        Vfmap{i}(voxel_spect{i}{j}.voxel(1),voxel_spect{i}{j}.voxel(2),voxel_spect{i}{j}.voxel(3),:) = Vf;
    end
end



%Vfmapnii = make_nii(Vfmap);
%save_nii(Vfmapnii,[datadir 'Vf_spect_maps_T1.nii.gz'])
%save([datadir 'Vf_spect_maps_T1'],'Vfmap');


%
for i=1:length(zslice)
    figure; hold on;
    for j=1:Ncomp
        subtightplot(1,Ncomp,j,0.02)
        imagesc(imrotate(Vfmap{i}(:,:,zslice(i),j),270))
        caxis([0 1])
        %title([num2str(T1_bounds(i)) ' < T1 < ' num2str(T1_bounds(i+1))])
        axis off
        %title(compartment_titles{j},'FontSize',20)
        colormap hot
    end
    set(gcf,'Position',[-106        1388        1499         397])
    
    if saveon
        print_to_formats(['/Users/paddyslator/Documents/OtherDocs/abstracts/Vfmap_' scans{i}],fig_formats)
    end
    
end

%plot a colorbar separately
figure;
imagesc(imrotate(Vfmap{i}(:,:,zslice(i),j),270))
caxis([0 1])
colormap hot    
colorbar('FontSize',20)

if saveon
   print_to_formats('/Users/paddyslator/Documents/OtherDocs/abstracts/Vfmap_colorbar',fig_formats)   
end


%% plot an example spectra from 

%find some voxels with high Vf for each component 
for i=1:length(scans)
    for j=1:Ncomp
        find(Vfmap{i}(:,:,:,j)>0.8)
        [vox1,vox2,vox3] = ind2sub(size(Vfmap{1}(:,:,:,1)),indices);
        l=1;
        for k=1:100
            if all(voxel_spect{i}{k}.voxel == [vox1 vox2 vox3])
                while l<5
                    figure;
                    plot_3D_spectrum(voxel_spect{i}{j}.F)
                    l=l+1;
                end
            end
        end
    end
end






%%
%fit spectrum in some other ROIs
maskfile = 'iron_mask.nii.gz';
maskfile = 'CSF_mask.nii.gz';

mask = load_nii(strcat(datadir, maskfile));
mask = logical(mask.img);

sig = calculate_mean_signal(input.img,mask);

output=ILT_3D(sig,gradechoinv,ILT_options);

plot_3D_spectrum(output)




%% 
datadir = '/Users/paddyslator/Dropbox/placentaJhu/Data_other/t1diff/zebra102/';

imgfile = '9N_r.nii.gz';

input = load_nii([datadir imgfile]);
input.img = double(input.img) + eps;
resortedTI = input.img;
dims = input.hdr.dime.dim(2:5);

gradechoinv = csvread([datadir 'gradechoinv_paddy_corrections.csv']);

maskfile = 'brain_mask.nii.gz';
mask = load_nii(strcat(datadir, maskfile));
mask = logical(mask.img);

Nx=dims(1); Ny=dims(2); Nz=dims(3);
Nd=dims(4); 






%% calculate the ROI signal

z=33;
slicemask = zeros(size(mask));
slicemask(:,:,z) = mask(:,:,z);
sigcor = calculate_mean_signal(inputcor.img,slicemask);
sig = calculate_mean_signal(input.img,slicemask);


%% fit the ROI spectrum
ILT_options.Nk1=20;
ILT_options.Nk2=20;
ILT_options.Nk3=20;

ILT_options.mink1 = 0.00005;
ILT_options.maxk1 = 0.01;
ILT_options.mink2 = 5;
ILT_options.maxk2 = 200;
ILT_options.mink3 = 250;
ILT_options.maxk3 = 5000;


ILT_options.alpha = 0.05;


output=ILT_3D(sig,gradechoinv,ILT_options);

% plot projections
figure;hold on;

i=1;
subtightplot(1,3,i,0.05,0.1, 0.05)
proj = squeeze(sum(output.F,1));
contour(output.w3,output.w2,proj)
xlabel('T1')
ylabel('T2*')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

i=2;
subtightplot(1,3,i,0.05,0.1, 0.05)
proj = squeeze(sum(output.F,2));
contour(output.w3,output.w1,proj)
xlabel('T1')
ylabel('D')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

i=3;
subtightplot(1,3,i,0.05,0.1, 0.05)
proj = squeeze(sum(output.F,3));
contour(output.w2,output.w1,proj)
xlabel('T2*')
ylabel('D')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')


