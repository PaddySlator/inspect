b1 = repmat([0 25 50 100 150 200 250 500 750 1000],[1 10]);
b2 = [zeros([1 10]) 25*ones(1,10) 50*ones(1,10) 100*ones(1,10) 150*ones(1,10)...
    200*ones(1,10) 250*ones(1,10) 500*ones(1,10) 750*ones(1,10) 1000*ones(1,10)];

gradechoinv = nan(length(b1),8);
gradechoinv(:,4) = b1;
gradechoinv(:,8) = b2;

%blood
d1 = 0.02;
%tissue
d2 = 0.002;

%do a super simple exchange simulation

%(could equally calculate these from exchange rates)
%fraction that goes from 1 to 2
f12 = 0.25;
%fraction that goes from 2 to 1
f21 = 0.25;
%fraction that stays in 1
f11 = 0.25;
%fraction that stays in 2
f22 = 1 - f11 - f21 - f12;

%signals for each of the populations
S12 = KernelDD([d1 d2], gradechoinv);
S21 = KernelDD([d2 d1], gradechoinv);
S11 = KernelDD([d1 d1], gradechoinv);
S22 = KernelDD([d2 d2], gradechoinv);

%total signal
S = f12*S12 + f21*S21 + f11*S11 + f22*S22;


kernelname = 'DD';
options = default_ILT_options(kernelname,gradechoinv);

%add a bit of regularisation
%options.reg = 1;
%options.alpha = 0.0001;

ILT_output = ILT(S,gradechoinv,options);

%plot the exchange spectrum
contour(ILT_output.grid{2},ILT_output.grid{1},ILT_output.F)



%now plot a big image with lots of stuff!
%image volume fractions
imgdim = [50 50 1];
%define mask
mask = ones(imgdim);

% S12map = zeros(imgdim);
% S21map = zeros(imgdim);
% S11map = zeros(imgdim);
% S22map = zeros(imgdim);

%set volume fractions across image
ncomp = 4;

vfimg = zeros([imgdim ncomp]);
vfimg(:,:,:,1) = repmat(linspace(0,0.5,imgdim(1)), [imgdim(2) 1])';
vfimg(:,:,:,2) = repmat(linspace(0,0.5,imgdim(1)),[imgdim(2) 1]);
vfimg(:,:,:,3) = abs(1 - sum(vfimg(:,:,:,1:2),4));
vfimg(:,:,:,4) = vfimg(:,:,:,3);
vfimg(:,:,:,3) = tril(vfimg(:,:,:,3),-1);
vfimg(:,:,:,4) = triu(vfimg(:,:,:,4));


%simulate a D-D experiment 

%load the gradechoinv file
gradechoinv_filename = 'inspect/examples/placenta_gradechoinv.txt';

%set the kernel option
kernel.name = 'DD';

%define the canonical spectral components
d1_comp = [d1 d2 d1 d2]; 
d2_comp = [d2 d1 d1 d2]; 
spectral_comp = [d1_comp;d2_comp];

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
            %add noise
            %S = add_noise(S,SNR,noisetype);
            %scale
            %S0=100;
            %S=S*S0;
            %assign to image
            simimg(x,y,z,:) = S;                 
        end
    end
end

%%% do the inspect and voxelwise fits %%%
%set the number of components to the ground truth
map_options.ncomp = 4;
%simple example so should converge quickly
map_options.maxiter = 5;
DD_siminspectmap = inspect_map(simimg,gradechoinv,mask,'DD',map_options);

%plot the ground truth maps
figure; hold on;
for i=1:4
    subplot(4,1,i);
    imagesc(vfimg(:,:,1,i))
end


%plot the ground truth components
figure; hold on;
plot(d1, d1,'x');
plot(d1, d2,'x');
plot(d2, d1,'x');
plot(d2, d2,'x');
xlim([DD_siminspectmap.options.ILT.mink(1) DD_siminspectmap.options.ILT.maxk(1)])
ylim([DD_siminspectmap.options.ILT.mink(2) DD_siminspectmap.options.ILT.maxk(2)])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')


%TO DO: fit inspect_vox (i.e. normal DEXSY) to compare







