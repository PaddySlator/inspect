
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

simoptions.SNR = 200;
simoptions.noisetype = 'rician';

ncomp = length(spectparams.T2);

vfimg = zeros([dim ncomp]);
simimg = zeros([dim size(gradechoinv,1)]);


% %set volume fractions across image
% vfimg(1:3,:,:,1) = 0.75;
% vfimg(1:3,:,:,2) = 0.25;
% vfimg(1:3,:,:,3) = 0;
%
% vfimg(4:8,:,:,1) = 0.5;
% vfimg(4:8,:,:,2) = 0.25;
% vfimg(4:8,:,:,3) = 0.25;
%
% %vfimg(7:8,:,:,1) = 0;
% %vfimg(7:8,:,:,2) = 1;
% %vfimg(7:8,:,:,3) = 0;
%
% vfimg(9:10,:,:,1) = 0.0;
% vfimg(9:10,:,:,2) = 0.25;
% vfimg(9:10,:,:,3) = 0.75;


%these are a bit more interesting!
%set volume fractions across image
clear vfimg
vfimg(:,:,:,1) = repmat(linspace(0,0.5,dim(1)), [dim(2) 1])';
vfimg(:,:,:,2) = repmat(linspace(0,0.5,dim(1)),[dim(2) 1]);
vfimg(:,:,:,3) = abs(1 - sum(vfimg(:,:,:,1:2),4));
vfimg(:,:,:,4) = vfimg(:,:,:,3);

vfimg(:,:,:,3) = tril(vfimg(:,:,:,3),-1);
vfimg(:,:,:,4) = triu(vfimg(:,:,:,4));






%even more interesting! loads of shapes, some with gradients, some solid blocks!
%circle

%square
% squarexpos = 5:15;
% squareypos = 5:15;
% vfimg(squarexpos,squareypos,:,1) = 0.5;
% vfimg(squarexpos,squareypos,:,2) = 0.25;
% vfimg(squarexpos,squareypos,:,3) = 0.5;
%
% squarexpos = 20:30;
% vfimg(squarexpos,squareypos,:,1) = 0.1;
% vfimg(squarexpos,squareypos,:,2) = 0.1;
% vfimg(squarexpos,squareypos,:,3) = 0.8;
%
% squarexpos = 35:45;
% vfimg(squarexpos,squareypos,:,1) = 0.1;
% vfimg(squarexpos,squareypos,:,2) = 0.8;
% vfimg(squarexpos,squareypos,:,3) = 0.1;

% = repmat(linspace(0,1,10), [1 5])';
% vfimg(:,:,:,2) = repmat(linspace(1,0,10), [1 5])';
% vfimg(:,:,:,3) = repmat(linspace(0,0.5,10), [1 5])';




%do some lines of different volume fractions
% comp1vfs = [0 0 0 0.25 0.25 0.25 0.25 0.5 0.5 0.5 0.75 0.75 ];
% comp2vfs = [0.25 0.5 0.75 0.25 0.25 0.75 0 0 0.5 0.25 0.25 0];
% comp3vfs = [0.75 0.5 0.25 0.5 0.5 0 0.75 0.5 0 0.25 0 0.25];
%
% compwidth = 5;
% for i=1:dim(1)/compwidth
%     vfimg(:, compwidth*(i-1)+1 : compwidth*i,:, 1 ) = comp1vfs(i);
%     vfimg(:, compwidth*(i-1)+1 : compwidth*i,:, 2 ) = comp2vfs(i);
%     vfimg(:, compwidth*(i-1)+1 : compwidth*i,:, 3 ) = comp3vfs(i);
% end


%normalise
vfimg = vfimg./sum(vfimg,4);



%randomly sample spectral volume fractions them all
% vfimg = rand(size(vfimg));
% %normalise
% for x=1:dim(1)
%     for y=1:dim(2)
%         for z=1:dim(3)
%             vfimg(x,y,z,:) = vfimg(x,y,z,:)./sum(vfimg(x,y,z,:)) ;
%         end
%     end
% end


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


denoise = 0;
if denoise
    nbins=10;
    [denoised_simimgvox, snoise, nsig] = MPmoments(simimgvox,nbins);

    %put denoised back into image form
    denoised_simimg = voxel_to_image(denoised_simimgvox,...
        voxind,...
        dim);

    rawsimimg = simimg;
    simimg = denoised_simimg;
else
    rawsimimg = simimg;
end

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




%inspect_options.ILT  = default_ILT_options('DT2');
%for the initial mean fit
%inspect_options.ILT_mean = default_ILT_options('DT2');

%inspect_options.n_comp = 4;

%inspect_options.maxiter = 2 ;
%inspect_options.init = 'meanspectrum';

%inspect_options.weightstol = 10^-3;

%inspect_options.onF = 1;
%inspect_options.onweights = 1;

%inspect_options.save = saveon;
%inspect_options.save_path = figuredir;
%inspect_options.dirname = dir_string;
%inspect_options.scan_names = {''};

%inspect_options.relabel = 0;

inspect_options.save=0;
inspect_options.maxiter=2;

rawsiminspectmap = inspect_map(rawsimimg,gradechoinv,mask,'DT2',inspect_options);

if denoise
    inspect_options.dirname = [dir_string(1:end-1) '_denoised/' ];
    siminspectmap = inspect_map(simimg,gradechoinv,mask,inspect_options);
end






%% fit the spectrum to each voxel individually
ILT_options.Nk1=50;
ILT_options.Nk2=50;

ILT_options.mink1 = 2 * 10^-4;
ILT_options.maxk1 = 5000 * 10^-3;
ILT_options.mink2 = 10 ;
ILT_options.maxk2 = 150;

ILT_options.alpha = 0;

Vfvoxelwise = zeros([dim ncomp]);

l=1;
for x=1:dim(1)
    for y=1:dim(2)
        for z=1:dim(3)
            sig = squeeze(simimg(x,y,z,:));
            %output.(change{1}){i}{j}{k}{x,y,z} = ILT_2D(sig, gradechoinv,ILT_options);

            %don't store this because it becomes massive
            output = ILT_2D(sig, gradechoinv,ILT_options);
            %calculate vf map

            spectparams.T2 = 10^3*[0.05 0.06 0.07 0.08] ;
            spectparams.D = [0.0002 0.003 0.05 0.2];



            boundadc = [0 0.001;0.001 0.01; 0.01 0.1; 0.1 10];
            boundt2 = 10^3 * [0 0.055;0.055 0.065;0.065 0.075; 0.075 0.2];

            %boundadc = [0 0.01;0.01 0.06; 0.06 10];
            %boundt2 = [0 0.065;0.065 0.075; 0.075 0.2];




            Vfvoxelwise(x,y,z,:) = integrate_spectrum_2D(output,boundadc,boundt2);

            disp(['completed voxel ' num2str(l) ' of ' num2str(prod(dim))])
            l=l+1;
        end

    end
end




%make voxelwise spectrum map into nifti
%Vfvoxelwisenii = make_nii(Vfvoxelwise);


if saveon
    niftiwrite(Vfvoxelwise, [figuredir dir_string 'Vfvoxelwise.nii.gz']);
    %save_nii(Vfvoxelwisenii, [figuredir dir_string 'Vfvoxelwise.nii.gz']);
end

%niftiwrite(Vfvoxelwise, [figuredir dir_string 'Vfvoxelwisematlab.nii.gz']);





%% plot the output

%get the grid to plot the spectra on
grid = getkernelgrid(rawsiminspectmap.options.ILT);
%get the output spectra
Fcomp = rawsiminspectmap.iter{end}{end}.Fcomp;
%get the output voxelwise spectral weights
imgweights = rawsiminspectmap.iter{end}{end}.imgweights;

for i=1:rawsiminspectmap.options.ncomp
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
subfigdim = [3 rawsiminspectmap.options.ncomp];
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

ILT_options = default_ILT_options('D');

ILT_options.reg = 0;
%regularisation parameter
ILT_options.alpha = 0;

%regularisation for the mean fit
inspect_options.ILT_mean = ILT_options;

%regularisation for the inspect map
inspect_options.ILT = ILT_options;

inspect_options.n_comp = 4;

inspect_options.maxiter = 2 ;
inspect_options.init = 'random';
inspect_options.init = 'kmeans';
inspect_options.init = 'meanspectrum';
%inspect_options.init = 'user';

inspect_options.parallel = 0;

inspect_options.onhill = 0;
inspect_options.fmincon = 1;
inspect_options.updateF = 0;

ILT_options.reg = 0;
%regularisation parameter
ILT_options.alpha = 0;


inspect_options.weightstol = 10^-3;


inspect_options.onF = 1;
inspect_options.onweights = 1;

inspect_options.save = saveon;
inspect_options.save_path = figuredir;
inspect_options.dirname = dir_string;
inspect_options.scan_names = {''};

inspect_options.relabel = 0;

D_siminspectmap = inspect_map(D_simimg,gradechoinv,mask,inspect_options);






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

%% fit the ILT voxelwise

ILT_options = default_ILT_options('DT2T1');

DT2T1_voxoutput = cell(dim);

l=1;
for x=1:dim(1)
    for y=1:dim(2)
        for z=1:dim(3)  
            
            sig = squeeze(DT2T1_simimg(x,y,z,:));
            %output.(change{1}){i}{j}{k}{x,y,z} = ILT_2D(sig, gradechoinv,ILT_options);

            tic;
            %don't store this because it becomes massive
            DT2T1_voxoutput{x,y,z} = ILT(sig, gradechoinv,ILT_options);
            time=toc;
            
            disp(['voxel ' num2str(l) ' of ' num2str(prod(dim)) ', it took ' num2str(time) ' seconds.'])
            l=l+1;
        end
    end
end
            

%%

disp('want to do some default settings!!')

ILT_options = default_ILT_options('DT2T1');

ILT_options.reg = 0;
%regularisation parameter
ILT_options.alpha = 0;


%regularisation for the mean fit
inspect_options.ILT_mean = ILT_options;


%regularisation for the inspect map
inspect_options.ILT = ILT_options;


inspect_options.n_comp = 4;

inspect_options.maxiter = 1 ;
inspect_options.init = 'meanspectrum';

inspect_options.fmincon = 1;

inspect_options.weightstol = 10^-3;

inspect_options.onF = 1;
inspect_options.onweights = 1;

inspect_options.save = saveon;

disp('pretty sure it relabels anyway!')
inspect_options.relabel = 0;


DT2T1_siminspectmap = inspect_map(DT2T1_simimg,gradechoinv,mask,inspect_options);



%% old stuff that might be useful

%% simulate using spectra instead

% %get rep spectra componenets
% spectparams.T2 = [0.03 0.08];
% spectparams.D = [0.002 0.02];
%
% spectparams.T2 = [0.06 0.10 0.14];
% spectparams.D = [0.003 0.02 0.2];
%
%
% spectparams.f = [1 0 0];
% S1 = simulate_multiexp_signal(spectparams,gradechoinv,SNR);
% spectparams.f = [0 1 0];
% S2 = simulate_multiexp_signal(spectparams,gradechoinv,SNR);
% spectparams.f = [0 0 1];
% S3 = simulate_multiexp_signal(spectparams,gradechoinv,SNR);
%
%
% output1 = ILT_2D(S1,gradechoinv,ILT_options);
% output2 = ILT_2D(S2,gradechoinv,ILT_options);
% output3 = ILT_2D(S3,gradechoinv,ILT_options);
%
% Fcomp{1} = output1.F;
% Fcomp{2} = output2.F;
% Fcomp{3} = output3.F;
%
% simimgspect = zeros([dim size(gradechoinv,1)]);
%
% l=1;
% for x=1:dim(1)
%     for y=1:dim(2)
%         for z=1:dim(3)
%             F = construct_spectrum_from_components(Fcomp, vfimg(x,y,z,:) );
%
%             S = simulate_spectrum_signal(F,gradechoinv,ILT_options,SNR);
%
%             simimgspect(x,y,z,:) = S;
%
%             disp(['voxel ' num2str(l) ' of ' num2str(prod(dim)) ' finished'])
%             l=l+1;
%         end
%     end
% end



%% test code that puts together spectral components
% %grid for calculating the spectrum on
% [W1,W2] = meshgrid(w1,w2);
%
% mu{1} = [w1(25) w2(25)];
% sigma{1} = [0.00001 0; 0 0.0001];
%
% Fcomp{1} = mvnpdf([W1(:) W2(:)], mu{1}, sigma{1});
% Fcomp{1} = reshape(Fcomp{1}, ILT_options.Nk2, ILT_options.Nk1);
%
%
% i=2;
% mu{i} = [w1(10) w2(10)];
% sigma{i} = [0.0000001 0; 0 0.0001];
% Fcomp{i} = mvnpdf([W1(:) W2(:)], mu{i}, sigma{i});
% Fcomp{i} = reshape(Fcomp{i}, ILT_options.Nk2, ILT_options.Nk1);
%
%
% for i=1:2
%     figure;hold on;
%     contour(w2,w1,Fcomp{i});
%     set(gca, 'XScale', 'log')
%     set(gca, 'YScale', 'log')
% end
%
%
%
% %put the components together
% F = construct_spectrum_from_components(Fcomp,[0.5 0.5]);
%
% figure; hold on;
% contour(w2,w1,F);
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
