clear

%set some options that are the same for all simulations

%saving options
saveon=0;
map_options.save=saveon;
vox_options.save=saveon;

%noise level and type 
SNR = 200;
noisetype = 'rician';


%% Simulate a simple volume fraction image.
% Each voxel is the weighting of the corresponding 
% canonical spectral component. All simulations use 
% the same volume fraction image

%image dimension
imgdim = [20 20 1];
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


%% Simulate a T2-D (equivalently T2*-D) experiment.
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
%preallocate some stuff
simimg = zeros([imgdim size(gradechoinv,1)]);

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

% T2-D example fitting


%
%%% do the inspect and voxelwise fits %%%
%set the number of components to the ground truth
map_options.ncomp = 4;
%simple example so should converge quickly
map_options.maxiter=1;
T2D_siminspectmap = inspect_map(simimg,gradechoinv,mask,'DT2',map_options);


% fit voxelwise spectra and integrate in spectral ROIs to get
% volume fraction maps

%set spectral ROIs
vox_options.sROI = [];
vox_options.sROI{1} = [0 0.001;0.001 0.01; 0.01 0.1; 0.1 10];
vox_options.sROI{2} = 10^3 * [0 0.055;0.055 0.065;0.065 0.075; 0.075 0.2];

%voxelwise fit
T2D_simvoxfit = inspect_vox(simimg,gradechoinv,mask,'DT2',vox_options);




%%% plot the output %%%

%unpack some output variables
%spectral grid
grid = T2D_siminspectmap.options.ILT.grid;
%output spectra
Fcomp = T2D_siminspectmap.iter{end}{end}.Fcomp;
%output voxelwise spectral weights
imgweights = T2D_siminspectmap.iter{end}{end}.imgweights;

%plot the spectral components 
for i=1:T2D_siminspectmap.options.ncomp
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
plot_map_vs_vox_sim(T2D_siminspectmap,T2D_simvoxfit,vfimg);






%% D example - e.g. model-free IVIM

%define an IVIM-type protocol
b = [0 50 100 200 300 400 600 800];
b=repmat(b,[1 3]);
%put it into gradechoinv format
gradechoinv = ones(length(b),4);
gradechoinv(:,4) = b;

%set the kernel option
kernel.name = 'D';

%simulate example spectrum
%diffusivities
d = [0.001 0.003 0.01 0.1];
spectral_comp = d;
       
%simulate the image
D_simimg = zeros([imgdim size(gradechoinv,1)]);%preallocate
for x=1:imgdim(1) %loops over voxels
    for y=1:imgdim(2)
        for z=1:imgdim(3)            
            %get volume fractions for this voxel
            f = vfimg(x,y,z,:);
            %simulate the signal by summing each component's signal              
            S = 0;            
            for i=1:ncomp
                %parameters for this component in this voxel
                kernel.params = spectral_comp(:,i);               
                S = S + f(i) * KernelMeas(kernel,gradechoinv);
            end                        
            %add noise
            S = add_noise(S,SNR,'rician');
            %assign
            D_simimg(x,y,z,:) = S;
        end
    end
end


%%% do the fits %%%
%set the number of components to the ground truth
map_options.ncomp = 4;
%simple example so should converge quickly
map_options.maxiter=1;
%continuous inspect 
D_siminspectmap = inspect_map(D_simimg,gradechoinv,mask,'D',map_options);

%voxelwise
%set spectral ROIs
vox_options.sROI = [];
vox_options.sROI{1} = [0 0.002;0.002 0.005; 0.005 0.05; 0.05 1];
%do the fit
D_simvoxfit = inspect_vox(D_simimg,gradechoinv,mask,'D',vox_options);



%%% plot the output %%%

%unpack some output variables
%spectral grid
grid = D_siminspectmap.options.ILT.grid;
%output spectra
Fcomp = D_siminspectmap.iter{end}{end}.Fcomp;
%output voxelwise spectral weights
imgweights = D_siminspectmap.iter{end}{end}.imgweights;

%plot the spectral components 
for i=1:D_siminspectmap.options.ncomp
    figure;hold on;
    plot(grid{1},Fcomp{i})
    
    plot(spectral_comp(:,i),0,'rx')
    
    title(['Component ' num2str(i)])
    legend({'InSpect fit','Ground Truth'})
    set(gca, 'XScale', 'log');
    xlabel('ADC (mm^2/s)')
end

%plot the maps
plot_map_vs_vox_sim(D_siminspectmap,D_simvoxfit,vfimg);




%% T1-T2-D example

%simulate ZEBRA brain type experiment (e.g. MUDI challenge
%http://cmic.cs.ucl.ac.uk/cdmri/challenge.html)

gradechoinv_filename = 'inspect/examples/mudi_gradechoinv.txt';
gradechoinv = load(gradechoinv_filename);


kernel.name = 'DT2T1';
%the canonical spectral components
d = [0.001 0.003 0.003 0.001];
t2 = [80 80 120 120];
t1 = [500 1000 1500 2000];


spectral_comp = [d;t2;t1];

%preallocate the simulated image
DT2T1_simimg = zeros([imgdim size(gradechoinv,1)]);
%simulate the signal over image voxels
for x=1:imgdim(1)
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
            %add noise
            S = add_noise(S,SNR,'rician');
            %assign to image
            DT2T1_simimg(x,y,z,:) = S;
        end
    end
end

% fit the ILT voxelwise and continuous inspect 

%set the number of components to the ground truth
map_options.ncomp = 4;
%simple example so should converge quickly
map_options.maxiter=1;
%DT2T1_siminspectmap = inspect_map(DT2T1_simimg,gradechoinv,mask,kernel.name,map_options);

%clear options
vox_options.sROI = [];
vox_options.sROI{1} = [0 0.002;0.002 0.01; 0.002 0.1; 0 0.002];
vox_options.sROI{2} = 10^3 * [0 0.1;0 0.1;0.1 0.2; 0.1 0.2];
vox_options.sROI{3} = [0 750;750 1250;1250 1750; 1750 5000];

DT2T1_simvoxfit = inspect_vox(DT2T1_simimg,gradechoinv,mask,kernel.name,vox_options);






% plot the output

%unpack some output variables
%spectral grid
grid = DT2T1_siminspectmap.options.ILT.grid;
%output spectra
Fcomp = DT2T1_siminspectmap.iter{end}{end}.Fcomp;
%output voxelwise spectral weights
imgweights = DT2T1_siminspectmap.iter{end}{end}.imgweights;

%plot the spectral components 
nproj = length(kernel.params);%get the number of projections/subplots
for i=1:DT2T1_siminspectmap.options.ncomp
    plot_multidim_spectrum(Fcomp{i},grid,kernel.name,spectral_comp(:,i))
      
    
    title(['Component ' num2str(i)])
    legend({'InSpect fit','Ground Truth'})
    set(gca, 'XScale', 'log');
    xlabel('ADC (mm^2/s)')
end

%plot the maps
plot_map_vs_vox_sim(DT2T1_siminspectmap,DT2T1_simvoxfit,vfimg);








