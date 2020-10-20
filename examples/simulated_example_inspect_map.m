
clear

addpath(genpath('inspect'))

%set some options that are the same for all simulations

%saving options
saveon=0;
inspect_options.save=saveon;

%noise level and type 
SNR = 400;
noisetype = 'rician';


%% Simulate volume fraction image.
% Each voxel is the weighting of the corresponding 
% canonical spectral component. All simulations use 
% the same volume fraction image

%image dimension
imgdim = [10 10 1];
%define mask
mask = ones(imgdim);
%number of spectral components 
ncomp = 4;

%set volume fractions across image
vfimg = zeros([imgdim ncomp]);
vfimg(:,:,:,1) = repmat(linspace(0,0.5,imgdim(1)), [imgdim(2) 1])';
vfimg(:,:,:,2) = repmat(linspace(0,0.5,imgdim(1)),[imgdim(2) 1]);
vfimg(:,:,:,3) = abs(1 - sum(vfimg(:,:,:,1:2),4));
vfimg(:,:,:,4) = vfimg(:,:,:,3);
vfimg(:,:,:,3) = tril(vfimg(:,:,:,3),-1);
vfimg(:,:,:,4) = triu(vfimg(:,:,:,4));

%normalise
vfimg = vfimg./sum(vfimg,4);


%% Simulate a T2-D (equiv. T2*-D) experiment.
% Canonical spectral component values and MR acquisition 
% parameters are based on those in: 
%
% Slator et al., Combined diffusion-relaxometry MRI to identify 
% dysfunction in the human placenta, MRM 2019.
% https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.27733
%

%load the gradechoinv file
gradechoinv_filename = 'inspect/examples/placenta_gradechoinv.txt';
gradechoinv = load(gradechoinv_filename);

%set the kernel option
kernel.name = 'DT2';

%define the canonical spectral components
d = [0.0002 0.003 0.05 0.2]; %ADC
t2 = [50 60 70 80]; %T2
spectral_comp = [d;t2];

%now simulate the image
simimg = zeros([imgdim size(gradechoinv,1)]);%preallocate
for x=1:imgdim(1) %loop over voxels
    for y=1:imgdim(2)
        for z=1:imgdim(3)            
            %get volume fractions for this voxel
            f = squeeze(vfimg(x,y,z,:)); 
            %simulate the signal by summing each component's signal   
            S=0;            
            for i=1:length(f) %loop over components                     
                %parameters for this component in this voxel
                kernel.params = spectral_comp(:,i);
                S = S + f(i) * KernelMeas(kernel,gradechoinv);                                                              
            end 
            %normalise
            b0teminindex = 1;           
            S=S./S(b0teminindex);
            %add noise
            S = add_noise(S,SNR,noisetype);
            %scale
            %S0=100;
            %S=S*S0;
            %assign to image
            simimg(x,y,z,:) = S;
        end
    end
end

simimgfixed = simimg;



%% simulate with continuous canonical spectral components

%define the canonical spectral components
meand = [0.0002 0.003 0.05 0.2]; %ADC
meant2 = [50 60 70 80]; %T2

meanspectral_comp = [meand;meant2];

vard = [0.000000002 0.0000003 0.00005 0.002];
vart2 = [5 5 5 5];



%spectral_comp = [d;t2];

ncomp = length(meand);
Fsim = cell(ncomp,1);

%define the grid for the components to be defined on
grid_options.mink = [2*10^-4  5];
grid_options.maxk = [5  200];
grid_options.Nk = [50 50];
grid_options.loggrid = [1 0];
grid_options.kernel='DT2';       
grid_options.reg=0;

grid = getkernelgrid(grid_options);
allgridcombs = ( combvec(grid{:}) )';

Fvec=cell(ncomp,1);
for i=1:ncomp
    Fvec{i} = mvnpdf(allgridcombs,meanspectral_comp(:,i)',[vard(i) 0; 0 vart2(i)]);
end
%normalise (on the discrete grid) so that spectra have same maximums
for i=1:ncomp
    Fvec{i} = Fvec{i}./sum(Fvec{i}(:));
end
%reshape vector to the correct dim
F=cell(ncomp,1);
for i=1:ncomp
    F{i} = reshape(Fvec{i}, grid_options.Nk);
end

figure; hold on
for i=1:ncomp
    contour(grid{2},grid{1},F{i})    
end
set(gca, 'YScale', 'log');


tic;
%now simulate the image
simimg = zeros([imgdim size(gradechoinv,1)]);%preallocate
for x=1:imgdim(1) %loop over voxels
    for y=1:imgdim(2)
        for z=1:imgdim(3)            
            %get volume fractions for this voxel
            f = squeeze(vfimg(x,y,z,:)); 
            %get the effective spectrum for this voxel
            Fvox = construct_spectrum_from_components(F,f);                        
            %normalise this spectrum
            
            
            %simulate the signal using this spectrum   
            S=simulate_spectrum_signal(Fvox,gradechoinv,grid_options); 
            
            %normalise
            b0teminindex = 1;           
            S=S./S(b0teminindex);
            %add noise
            S = add_noise(S,SNR,noisetype);
            %scale
            %S0=100;
            %S=S*S0;
            %assign to image
            simimg(x,y,z,:) = S;
        end
    end
end

%S = simulate_spectrum_signal(F,gradechoinv,options,SNR)

toc;





%% set up everything for saving
if saveon
    %where to save the figures
    %figuredir = '/Users/paddyslator/Dropbox/PlacentaDocs/papers/inspect_map/miccai2020/simulations/';
    %figuredir = '/home/pslator/IPMI2019/simulations/';
    figuredir = pwd;
    
    %make a nice string for the directory
    dir_string = ['SNR_' num2str(SNR)];
    
    dir_string = [dir_string '_T2'];
    for i=1:length(t2)
        dir_string = [dir_string '_' num2str(t2(i))];
    end
    
    dir_string = [dir_string '_D'];
    for i=1:length(d)
        dir_string = [dir_string '_' num2str(d(i))];
    end
    
    dir_string = [dir_string '/'];
    
    mkdir([figuredir dir_string]);
end

% save the ground truth as a nifti
if saveon
    niftiwrite(vfimg, [figuredir dir_string 'Vf_ground_truth.nii.gz'],'Compressed',true);
end

%save the simulated image as nifti
if saveon
    niftiwrite(simimg, [figuredir dir_string 'simimg'],'Compressed',true);
end

if saveon
    inspect_options.save_path = figuredir;
    inspect_options.dirname = dir_string;
end

%%
%%% do the inspect and voxelwise fits %%%




%choose these options (default is pretty similar but these line up the plots
%nicely)
inspect_options.ILT.mink = [2*10^-4  5];
inspect_options.ILT.maxk = [5  200];

inspect_options.ILT_mean.mink = [2*10^-4  5];
inspect_options.ILT_mean.maxk = [5  200];

inspect_options.maxiter = 5;

inspect_options.sumto1 = 0;

% fit inspect continuous version
comps = 4;
inspect_options.ncomp = comps;
     
siminspectmap = inspect_map(simimg,gradechoinv,mask,'DT2',inspect_options);

% fit voxelwise spectra and integrate in spectral ROIs to get
% volume fraction maps

%set spectral ROIs
clear options
vox_options.sROI{1} = [0 0.001;0.001 0.01; 0.01 0.1; 0.1 10];
vox_options.sROI{2} = 10^3 * [0 0.055;0.055 0.065;0.065 0.075; 0.075 0.2];

vox_options.ILT = inspect_options.ILT;
vox_options.ILT.alpha=0.01;

vox_options.save=saveon;

if saveon
    vox_options.save_path = inspect_options.save_path;
    vox_options.dirname = inspect_options.dirname;
end

%voxelwise fit
simvoxfit = inspect_vox(simimg,gradechoinv,mask,'DT2',vox_options);


%%% plot the output %%%

%unpack some output variables
%spectral grid
grid = siminspectmap.options.ILT.grid;
%output spectra
Fcomp = siminspectmap.iter{end}{end}.Fcomp;
%output voxelwise spectral weights
imgweights = siminspectmap.iter{end}{end}.imgweights;

%plot the spectral components 
for i=1:siminspectmap.options.ncomp
    figure;hold on;
    contour(grid{2},grid{1},Fcomp{i})
    
    plot(spectral_comp(2,i),spectral_comp(1,i),'rx')
    title(['Component ' num2str(i)])
    legend({'InSpect fit','Ground Truth'})
    set(gca, 'YScale', 'log');
    xlabel('T2* (s)')
    ylabel('ADC (mm^2/s)')
end

%plot a comparison of the volume fraction maps
plot_map_vs_vox_sim(siminspectmap,simvoxfit,vfimg);



return
   



%% D example - i.e. model-free IVIM

% %define an IVIM-type protocol
% b = [0 50 100 200 300 400 600 800];
% b=repmat(b,[1 3]);
% %put it into gradechoinv format
% gradechoinv = ones(length(b),4);
% gradechoinv(:,4) = b;
% 
% %set the kernel option
% kernel.name = 'D';
% 
% %simulate example spectrum
% %diffusivities
% d = [0.001 0.003 0.01 0.1];
% spectral_comp = d;
%        
% %simulate the image
% D_simimg = zeros([imgdim size(gradechoinv,1)]);%preallocate
% for x=1:imgdim(1) %loops over voxels
%     for y=1:imgdim(2)
%         for z=1:imgdim(3)            
%             %get volume fractions for this voxel
%             f = vfimg(x,y,z,:);
%             %simulate the signal by summing each component's signal              
%             S = 0;            
%             for i=1:ncomp
%                 %parameters for this component in this voxel
%                 kernel.params = spectral_comp(:,i);               
%                 S = S + f(i) * KernelMeas(kernel,gradechoinv);
%             end                        
%             %add noise
%             S = add_noise(S,SNR,'rician');
%             %assign
%             D_simimg(x,y,z,:) = S;
%         end
%     end
% end
% 
% 
% %%% do the fits %%%
% 
% %continuous inspect 
% D_siminspectmap = inspect_map(D_simimg,gradechoinv,mask,'D',inspect_options);
% 
% %voxelwise
% %set spectral ROIs
% clear options
% options.sROI{1} = [0 0.002;0.002 0.005; 0.005 0.05; 0.05 1];
% %do the fit
% D_simvoxfit = inspect_vox(D_simimg,gradechoinv,mask,'D',options);
% 
% 
% 
% %%% plot the output %%%
% 
% %unpack some output variables
% %spectral grid
% grid = D_siminspectmap.options.ILT.grid;
% %output spectra
% Fcomp = D_siminspectmap.iter{end}{end}.Fcomp;
% %output voxelwise spectral weights
% imgweights = D_siminspectmap.iter{end}{end}.imgweights;
% 
% %plot the spectral components 
% for i=1:D_siminspectmap.options.ncomp
%     figure;hold on;
%     plot(grid{1},Fcomp{i})
%     
%     plot(spectral_comp(:,i),0,'rx')
%     
%     title(['Component ' num2str(i)])
%     legend({'InSpect fit','Ground Truth'})
%     set(gca, 'XScale', 'log');
%     xlabel('ADC (mm^2/s)')
% end
% 
% %plot the maps
% plot_map_vs_vox_sim(D_siminspectmap,D_simvoxfit,vfimg);


%% T1-T2-D example

% simulate ZEBRA brain type experiment (e.g. MUDI challenge
% http://cmic.cs.ucl.ac.uk/cdmri/challenge.html)

% gradechoinv_filename = 'inspect/examples/mudi_gradechoinv.txt';
% gradechoinv = load(gradechoinv_filename);
% 
% 
% kernel.name = 'DT2T1';
% %the canonical spectral components
% d = [0.002 0.001 0.003 0.002];
% t2 = [80 80 120 120];
% t1 = [500 1000 1500 2000];
% 
% spectral_comp = [d;t2;t1];
% 
% %preallocate the simulated image
% DT2T1_simimg = zeros([imgdim size(gradechoinv,1)]);
% %simulate the signal over image voxels
% for x=1:imgdim(1)
%     for y=1:imgdim(2)
%         for z=1:imgdim(3)            
%             %get volume fractions for this voxel
%             f = squeeze(vfimg(x,y,z,:));                        
%             %simulate the signal by summing each component's signal   
%             S=0;
%             for i=1:length(f) %loop over components                     
%                 %parameters for this component in this voxel
%                 kernel.params = spectral_comp(:,i);
%                 S = S + f(i) * KernelMeas(kernel,gradechoinv);                                                              
%             end     
%             %add noise
%             S = add_noise(S,SNR,'rician');
%             %assign to image
%             DT2T1_simimg(x,y,z,:) = S;
%         end
%     end
% end
% 
% % fit the ILT voxelwise and continuous inspect 
% 
% %clear options
% clear options
% options.sROI{1} = [0 0.001;0.001 0.01; 0.01 0.1; 0.1 10];
% options.sROI{2} = 10^3 * [0 0.055;0.055 0.065;0.065 0.075; 0.075 0.2];
% options.sROI{3} = [0 500;500 1000;1000 2000; 2000 5000];
% 
% DT2T1_simvoxfit = inspect_vox(DT2T1_simimg,gradechoinv,mask,kernel.name);
% 
% DT2T1_siminspectmap = inspect_map(DT2T1_simimg,gradechoinv,mask,kernel.name,inspect_options);
% 
% 
% 
% % plot the output
% 
% %unpack some output variables
% %spectral grid
% grid = DT2T1_siminspectmap.options.ILT.grid;
% %output spectra
% Fcomp = DT2T1_siminspectmap.iter{end}{end}.Fcomp;
% %output voxelwise spectral weights
% imgweights = DT2T1_siminspectmap.iter{end}{end}.imgweights;
% 
% %plot the spectral components 
% nproj = length(kernel.params);%get the number of projections/subplots
% for i=1:DT2T1_siminspectmap.options.ncomp
%     plot_multidim_spectrum(Fcomp{i},grid,kernel.name,spectral_comp(:,i))
%       
%     
%     title(['Component ' num2str(i)])
%     legend({'InSpect fit','Ground Truth'})
%     set(gca, 'XScale', 'log');
%     xlabel('ADC (mm^2/s)')
% end
% 
% %plot the maps
% plot_map_vs_vox_sim(DT2T1_siminspectmap,DT2T1_simvoxfit,vfimg);






