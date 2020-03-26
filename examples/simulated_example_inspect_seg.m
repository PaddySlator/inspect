




%% T2-D example



%simulate a simple image where each voxel has an associated spectrum 
dim = [10 10 10];


%load a placenta gradechoinv file
gradechoinv = load('/Users/paddyslator/Dropbox/t2sdiff/pip0111/grad_echo_inv.txt');
%want this in ms
gradechoinv(:,5) = 10^3 * gradechoinv(:,5);

nmeas = size(gradechoinv,1);



%associated spectrum values for each ROI
%[T2 D f]
% spectparams = cell(nroi,1);
% 
% spectparams{1}.T2 = [0.03 0.03 0.03];
% spectparams{1}.D = [0.002 0.002 0.002];
% spectparams{1}.f = [0.4 0.2 0.4];
% 
% spectparams{2}.T2 = [0.04 0.05 0.06];
% spectparams{2}.D = [0.002 0.03 0.3];
% spectparams{2}.f = [0.4 0.2 0.4];
% 
% spectparams{3}.T2 = [0.08 0.08 0.08];
% spectparams{3}.D = [0.002 0.03 0.3];
% spectparams{3}.f = [0.4 0.2 0.4];
% 
% 
% spectparams{1}.T2 = [0.03 0.04];
% spectparams{1}.D = [0.001 0.01];
% spectparams{1}.f = [0.7 0.3];
% 
% spectparams{2}.T2 = [0.03 0.04 0.06];
% spectparams{2}.D = [0.001 0.01 0.1];
% spectparams{2}.f = [0.5 0.25 0.25];
% 
% spectparams{3}.T2 = [0.03 0.04 0.06];
% spectparams{3}.D = [0.001 0.01 0.03];
% spectparams{3}.f = [0.2 0.4 0.4];


clear spectparams params
spectparams.T2 = 10^3 * [0.05 0.1 0.05 0.1];
spectparams.D = [0.001 0.001 0.02 0.02];
spectparams.f =  [1 1 1 1];


nroi = 4;

SNR=100;

simoptions.SNR = SNR;
simoptions.noisetype = 'rician';


%make the segmented image
simroiimg = zeros(dim);

simroiimg(1:2,:,:) = 1;
simroiimg(3:4,:,:) = 2;
simroiimg(5:7,:,:) = 3;
simroiimg(8:10,:,:) = 4;


%now simulate the image
simimg = zeros([dim nmeas]);
simweightsimg = zeros(dim);

for x=1:dim(1)
    for y=1:dim(2)
        for z=1:dim(3)
            %get the ROI
            roi = simroiimg(x,y,z);
            
            params.D = spectparams.D(roi);
            params.T2 = spectparams.T2(roi);
            params.f = spectparams.f(roi);
            
            S = simulate_multiexp_signal(params,gradechoinv,simoptions);
            
            simimg(x,y,z,:) = S;
            simweightsimg(x,y,z,roi) = 1 - eps;
        end
    end
end




%% fit inspect segmentation model
mask = ones(dim);

%inspect segmentation version

% ILT_options = default_ILT_options('DT2');
% 
% inspect_seg_options.ILT = ILT_options;
% 
% inspect_seg_options.ILT.alpha=0.001;
% 
% 
% %the number of clusters is known from the simulation
% inspect_seg_options.EM.n_clusters = nroi;
% inspect_seg_options.EM.n_steps = 10;
% inspect_seg_options.EM.init = 'kmeans';
% 
% inspect_seg_options.ILT.kernel = 'DT2';
% 
% inspect_seg_options.save = 1;


t2d_output = inspect_seg(simimg,gradechoinv,mask,'DT2');


return 




%% D-D example






%% D example

%simulate an IVIM image

saveon=0;

%simulate a simple image where each voxel has an associated spectrum 
dim = [10 10 1];

%define the protocol
b = [0 50 100 300 600];
%b = [0 10 20 30 40 50 100 150 200 300 400 600 800];
%b = 0:800;

b=repmat(b,[1 3]);

gradechoinv = ones(length(b),5);
gradechoinv(:,4) = b;
%TE
gradechoinv(:,5) = 74;

SNR=200;

%number of ROIs - each with a different associated spectrum
nroi = 3;

%
D = [0.001 0.001 0.001];
Dv = [0.1 0.1 0.01];
f = [0.1 0.2 0.3];

S = zeros(length(b),1);

simroiimg = zeros(dim);

simroiimg(1:round(dim(1)/3),:,:) = 1;
simroiimg((1+round(dim(1)/3)):round(2 * dim(1)/3),:,:) = 2;
simroiimg((1+round(2 * dim(1)/3)):dim(1),:,:) = 3;


simimg = zeros([dim length(b)]);

    
for x=1:dim(1)
    for y=1:dim(2)
        for z=1:dim(3)
            roi = simroiimg(x,y,z);
            
            for i=1:length(b)
                S(i) = f(roi) * exp(-b(i) * Dv(roi)) ...
                    + (1 - f(roi)) * exp(-b(i) * D(roi));  
                                
            end
            
            S = add_noise(S,SNR,'rician');
                                         
            simimg(x,y,z,:) = S;
        end
    end
end


%% fit the spectrum ILT voxelwise

ILT_options.Nk = 200;
ILT_options.mink = 2*10^-4;
ILT_options.maxk = 5;
ILT_options.reg = 1;
ILT_options.alpha = 0.01;
ILT_options.kernel = 'D';
ILT_options.loggrid = 1;

voxelwise_output = cell(dim);

Vfvoxelwise = zeros([dim nroi]);

l=1;
for x=1:dim(1)
    for y=1:dim(2)
        for z=1:dim(3)
            sig = squeeze(simimg(x,y,z,:));
            %output.(change{1}){i}{j}{k}{x,y,z} = ILT_2D(sig, gradechoinv,ILT_options);

            %don't store this because it becomes massive
            voxelwise_output{x,y,z} = ILT(sig, gradechoinv,ILT_options);
            %calculate vf map
            spectparams.D = [0.0002 0.003 0.05 0.2];

            boundadc = [0 0.001;0.001 0.01; 0.01 0.1; 0.1 10];

            %Vfvoxelwise(x,y,z,:) = integrate_spectrum_2D(voxelwise_output{x,y,z},boundadc,boundt2);

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





%% fit inspect segmentation model
mask = ones(dim);

%inspect segmentation version
% ILT_options.Nk = 200;
% ILT_options.mink = 0;
% ILT_options.maxk = 0.1;
% ILT_options.alpha = 0.01;
% inspect_options.ILT.kernel = 'D';
%inspect_options.ILT.reg = 1;
%inspect_options.ILT.loggrid = 1;

ILT_options = default_ILT_options('D');

inspect_seg_options.ILT = ILT_options;

%the number of clusters is known from the simulation
inspect_seg_options.EM.n_clusters = nroi;
inspect_seg_options.EM.n_steps = 10;
inspect_seg_options.EM.init = 'kmeans';
inspect_seg_options.save = 0;
%inspect_options.save_path = [paperpath '/simulations/'];


d_output = inspect_seg(simimg,gradechoinv,mask,'D',inspect_seg_options);





%% T1-T2-D example

dim = [10 10 1];

gradechoinv = load('~/Dropbox/challenge/parameters.txt');
gradechoinv(:,7) = 7500;
%these are the other way around!
gradechoinvtemp = gradechoinv;
gradechoinv(:,5) = gradechoinvtemp(:,6);
gradechoinv(:,6) = gradechoinvtemp(:,5);

d = [0.002 0.001 0.003 0.002];
k = [0.7 1 0 0];
t2 = [100 100 150 150];
t1 = [500 1000 1500 2000];

SNR = 10000;

nroi = 4;

simroiimg = ones(dim);

simroiimg(1:2,:,:) = 1;
simroiimg(3:4,:,:) = 2;
simroiimg(5:7,:,:) = 3;
simroiimg(8:10,:,:) = 4;


simimg = zeros([dim size(gradechoinv,1)]);

for x=1:dim(1)
    for y=1:dim(2)
        for z=1:dim(3)
            roi = simroiimg(x,y,z);
            
            %params = [d(roi) t2(roi) t1(roi)];
            params = [d(roi) k(roi) t2(roi) t1(roi)];
            
            %S = KernelDT2T1(params,gradechoinv);
            S = KernelDKT2T1(params,gradechoinv);
    
            S = add_noise(S,SNR,'rician');
                                         
            simimg(x,y,z,:) = S;
        end
    end
end


%% fit inspect segmentation model
mask = ones(dim);

%inspect segmentation version
% ILT_options.Nk = 200;
% ILT_options.mink = 0;
% ILT_options.maxk = 0.1;
% ILT_options.alpha = 0.01;
% inspect_options.ILT.kernel = 'D';
%inspect_options.ILT.reg = 1;
%inspect_options.ILT.loggrid = 1;

clear inspect_seg_options

%ILT_options = default_ILT_options('DKT2T1');
%inspect_seg_options.ILT = ILT_options;



%the number of clusters is known from the simulation
inspect_seg_options.nclus = 4;
inspect_seg_options.nstep = 5;

inspect_seg_options.maxiter = 10;

%inspect_seg_options.EM.init = 'kmeans';
%inspect_seg_options.save = 0;
%inspect_options.save_path = [paperpath '/simulations/'];


dkt2t1_output = inspect_seg(simimg,gradechoinv,mask,'DKT2T1',inspect_seg_options);






