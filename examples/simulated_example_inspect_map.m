
saveon=0;

%simulate a simple image where each voxel has an associated spectrum
dim = [10 10 1];


%load a placenta gradechoinv file
gradechoinv = load('/Users/paddyslator/Dropbox/t2sdiff/pip0111/grad_echo_inv.txt');
%want this in ms
gradechoinv(:,5) = gradechoinv(:,5)*10^3;

nmeas = size(gradechoinv,1);


%make some spectral components
ILT_options.Nk1=50;
ILT_options.Nk2=50;

ILT_options.Nk = [50 50];

ILT_options.mink1 = 2 * 10^-4;
ILT_options.maxk1 = 5000 * 10^-3;
ILT_options.mink2 = 10 * 10^-3;
ILT_options.maxk2 = 150 * 10^-3;

ILT_options.mink = [2*10^-4  10*10^-3];
%ILT_options.maxk = [5000*10^-3  150*10^-3];
ILT_options.maxk = [5000*10^-3  150*10^-3];

ILT_options.reg = 0;
ILT_options.alpha = 0;


w1 = logspace(log10(ILT_options.mink1),log10(ILT_options.maxk1), ILT_options.Nk1);
w2 = linspace(ILT_options.mink2,ILT_options.maxk2, ILT_options.Nk2);



%simulate example spectrum
%the T2 and D are the same, but each voxel has different volume fractions
%of these components - the "canonical spectrum components" are
spectparams.T2 = [0.06 0.07 0.08] ;

spectparams.D = [0.002 0.02 0.1];


spectparams.T2 = 10^3 * [0.05 0.06 0.07 0.08] ;

spectparams.D = [0.0002 0.003 0.05 0.2];



%spectparams.T2 = [0.03 0.08];
%spectparams.D = [0.002 0.02];

simoptions.SNR = 25;
simoptions.noisetype = 'rician';

ncomp = length(spectparams.T2);

vfimg = zeros([dim ncomp]);
simimg = zeros([dim size(gradechoinv,1)]);



%these are a bit more interesting!
%set volume fractions across image
clear vfimg
vfimg(:,:,:,1) = repmat(linspace(0,0.5,dim(1)), [dim(2) 1])';
vfimg(:,:,:,2) = repmat(linspace(0,0.5,dim(1)),[dim(2) 1]);
vfimg(:,:,:,3) = abs(1 - sum(vfimg(:,:,:,1:2),4));
vfimg(:,:,:,4) = vfimg(:,:,:,3);

vfimg(:,:,:,3) = tril(vfimg(:,:,:,3),-1);
vfimg(:,:,:,4) = triu(vfimg(:,:,:,4));



%normalise
vfimg = vfimg./sum(vfimg,4);


for x=1:dim(1)
    for y=1:dim(2)
        for z=1:dim(3)
            %simulate the volume fractions for this voxel
%             %weightings for the vfs
%             vfweights = [1 2 3]';
%             f = vfweights .* rand(ncomp,1);
%             %normalise so they sum to 1
%             f = f./sum(f);
%             vfimg(x,y,z,:) = f;

            %get volume fractions for this voxel
            spectparams.f = vfimg(x,y,z,:);

            S0 = 1;

            S = S0.*simulate_multiexp_signal(spectparams,gradechoinv,simoptions);

            simimg(x,y,z,:) = S;
        end
    end
end




%put image into voxel form
[simimgvox, voxind] = image_to_voxel(simimg);


rawsimimg = simimg;

%% set up everything for saving

%where to save the figures
figuredir = '/Users/paddyslator/Dropbox/PlacentaDocs/papers/inspect_map/simulations/';

%make a nice string for the directory
dir_string = ['SNR_' num2str(simoptions.SNR)];

dir_string = [dir_string '_T2'];
for i=1:length(spectparams.T2)
    dir_string = [dir_string '_' num2str(spectparams.T2(i))];
end

dir_string = [dir_string '_D'];
for i=1:length(spectparams.D)
    dir_string = [dir_string '_' num2str(spectparams.D(i))];
end

dir_string = [dir_string '/'];

mkdir([figuredir dir_string]);



%% save the ground truth as a nifti

%make voxelwise spectrum map into nifti
%Vfimgnii = make_nii(vfimg);

if saveon
    niftiwrite(vfimg, [figuredir dir_string 'Vf_ground_truth.nii.gz'],'Compressed',true);
    %save_nii(Vfimgnii, [figuredir dir_string 'Vf_ground_truth.nii.gz']);
end



%save the simulated image as nifti
if saveon
    niftiwrite(simimg, [figuredir dir_string 'simimg'],'Compressed',true);
end



%% inspect continuous version
mask = ones(dim);

inspect_options.save=0;
inspect_options.maxiter=2;

siminspectmap = inspect_map(simimg,gradechoinv,mask,'DT2',inspect_options);



%% fit the spectrum to each voxel individually

options.sROI{1} = [0 0.001;0.001 0.01; 0.01 0.1; 0.1 10];
options.sROI{2} = 10^3 * [0 0.055;0.055 0.065;0.065 0.075; 0.075 0.2];

simvoxfit = fit_vox_spectra(simimg,gradechoinv,mask,'DT2',options);


   




%% plot the output

%get the grid to plot the spectra on
grid = getkernelgrid(siminspectmap.options.ILT);
%get the output spectra
Fcomp = siminspectmap.iter{end}{end}.Fcomp;
%get the output voxelwise spectral weights
imgweights = siminspectmap.iter{end}{end}.imgweights;

for i=1:siminspectmap.options.ncomp
    figure;hold on;
    contour(grid{2},grid{1},Fcomp{i})
    
    plot(spectparams.T2(i),spectparams.D(i),'rx')
    title(['Component ' num2str(i)])
    legend({'InSpect fit','Ground Truth'})
    set(gca, 'YScale', 'log');
    xlabel('T2* (s)')
    ylabel('ADC (mm^2/s)')
end

%
figure; hold on;
subfigdim = [3 siminspectmap.options.ncomp];
maporder = 1:subfigdim(2);
%cmax = [1 0.5 0.5];

for i=1:ncomp
    %subtightplot(subfigdim(1),subfigdim(2),i)
    subplot(subfigdim(1),subfigdim(2),i)
    imagesc(imgweights{1}(:,:,1,maporder(i)))
    
    colorbar
    %caxis([0 cmax(i)])
    if i==1
        ylabel('InSpect Maps')
    end

    %subtightplot(subfigdim(1),subfigdim(2),i+ncomp)
    subplot(subfigdim(1),subfigdim(2),i + ncomp)
    imagesc(vfimg(:,:,1,i))
    colorbar;
    if i==1
        ylabel('Ground truth volume fraction')
    end

    %subtightplot(subfigdim(1),subfigdim(2),i + 2 * ncomp)
    subplot(subfigdim(1),subfigdim(2),i + 2 * ncomp)
    imagesc(Vfvoxelwise(:,:,:,maporder(i)))
    colorbar
    %caxis([0 cmax(i)])

    if i==1
        ylabel('Voxelwise Maps')
    end
    
    
    
end

set(gcf,'Position',[-56 1011 1261 677])


return




%% D example

%simulate an IVIM image

%define the protocol
b = [0 50 100 300 600];

b=0:3000;
b=repmat(b,[1 3]);

gradechoinv = ones(length(b),5);
gradechoinv(:,4) = b;
%TE
gradechoinv(:,5) = 74;

SNR=2000;

%simulate example spectrum
%diffusivities
D = [0.001 0.003 0.01 0.1];

%simulate using the same volume fraction image
D_simimg = zeros([dim size(gradechoinv,1)]);

for x=1:dim(1)
    for y=1:dim(2)
        for z=1:dim(3)
            %simulate the volume fractions for this voxel
%             %weightings for the vfs
%             vfweights = [1 2 3]';
%             f = vfweights .* rand(ncomp,1);
%             %normalise so they sum to 1
%             f = f./sum(f);
%             vfimg(x,y,z,:) = f;

            %get volume fractions for this voxel
            f = vfimg(x,y,z,:);

            S0 = 1;
            
            S = 0;
            
            for j=1:ncomp
                S = S + f(j) * exp(-b * D(j));                  
            end
            
            S = S0 .* S;
            
            S = add_noise(S,SNR,'rician');

            D_simimg(x,y,z,:) = S;
        end
    end
end




%% fit inspect mapping version

% ILT_options = default_ILT_options('D');
% 
% ILT_options.reg = 0;
% %regularisation parameter
% ILT_options.alpha = 0;
% 
% %regularisation for the mean fit
% inspect_options.ILT_mean = ILT_options;
% 
% %regularisation for the inspect map
% inspect_options.ILT = ILT_options;
% 
% inspect_options.n_comp = 4;
% 
% inspect_options.maxiter = 2 ;
% inspect_options.init = 'random';
% inspect_options.init = 'kmeans';
% inspect_options.init = 'meanspectrum';
% %inspect_options.init = 'user';
% 
% inspect_options.parallel = 0;
% 
% inspect_options.onhill = 0;
% inspect_options.fmincon = 1;
% inspect_options.updateF = 0;
% 
% ILT_options.reg = 0;
% %regularisation parameter
% ILT_options.alpha = 0;
% 
% 
% inspect_options.weightstol = 10^-3;
% 
% 
% inspect_options.onF = 1;
% inspect_options.onweights = 1;
% 
% inspect_options.save = saveon;
% inspect_options.save_path = figuredir;
% inspect_options.dirname = dir_string;
% inspect_options.scan_names = {''};
% 
% inspect_options.relabel = 0;

D_siminspectmap = inspect_map(D_simimg,gradechoinv,mask,'D');

%fit voxelwise
D_voxfit = fit_vox_spectra(D_simimg,gradechoinv,mask,'D');







%% T1-T2-D example

dim = [10 10 1];

gradechoinv = load('~/Dropbox/challenge/parameters_new.txt');
gradechoinv(:,7) = 7500;
%these are the other way around!
gradechoinvtemp = gradechoinv;
gradechoinv(:,5) = gradechoinvtemp(:,6);
gradechoinv(:,6) = gradechoinvtemp(:,5);

d = [0.002 0.001 0.003 0.002];
t2 = [100 100 150 150];
t1 = [500 1000 1500 2000];

SNR = 10000;


%simulate using the same volume fraction image
DT2T1_simimg = zeros([dim size(gradechoinv,1)]);


for x=1:dim(1)
    for y=1:dim(2)
        for z=1:dim(3)            
            %get volume fractions for this voxel
            f = squeeze(vfimg(x,y,z,:));                        
               
            S=0;
            for i=1:length(f)
                S = S + f(i)*(KernelD(d(i),gradechoinv)...
                    .*KernelT2(t2(i),gradechoinv)...
                    .*KernelT1inv(t1(i),gradechoinv));
            end            
    
            S = add_noise(S,SNR,'rician');
                                         
            DT2T1_simimg(x,y,z,:) = S;
        end
    end
end

%% fit the ILT voxelwise and continuous inspect 


%options.nSROI = 4;
%options.sROI{1} = [0 0.001;0.001 0.01; 0.01 0.1; 0.1 10];
%options.sROI{2} = 10^3 * [0 0.055;0.055 0.065;0.065 0.075; 0.075 0.2];
%options.sROI{3} = [0 500;500 1000;1000 2000; 2000 5000];

simvoxfit = fit_vox_spectra(DT2T1_simimg,gradechoinv,mask,'DT2T1');

DT2T1_siminspectmap = inspect_map(DT2T1_simimg,gradechoinv,mask,'DT2T1');











