
%% T2-D example



%simulate a simple image where each voxel has an associated spectrum 
dim = [10 10 10];


%load a placenta gradechoinv file
gradechoinv = load('inspect/examples/placenta_gradechoinv.txt');

% %want this in ms
% gradechoinv(:,5) = 10^3 * gradechoinv(:,5);

nmeas = size(gradechoinv,1);



%associated spectrum values for each ROI
%[T2 D f]
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




% fit inspect segmentation model
mask = ones(dim);

t2d_output = inspect_seg(simimg,gradechoinv,mask,'DT2');

plot_inspect_seg(t2d_output)


return 










% %% D example
% 
% %simulate an IVIM image
% 
% saveon=0;
% 
% %simulate a simple image where each voxel has an associated spectrum 
% dim = [10 10 1];
% 
% %define the protocol
% b = [0 50 100 300 600];
% %b = [0 10 20 30 40 50 100 150 200 300 400 600 800];
% %b = 0:800;
% 
% b=repmat(b,[1 3]);
% 
% gradechoinv = ones(length(b),5);
% gradechoinv(:,4) = b;
% %TE
% gradechoinv(:,5) = 74;
% 
% SNR=200;
% 
% %number of ROIs - each with a different associated spectrum
% nroi = 3;
% 
% %
% D = [0.001 0.001 0.001];
% Dv = [0.1 0.1 0.01];
% f = [0.1 0.2 0.3];
% 
% S = zeros(length(b),1);
% 
% simroiimg = zeros(dim);
% 
% simroiimg(1:round(dim(1)/3),:,:) = 1;
% simroiimg((1+round(dim(1)/3)):round(2 * dim(1)/3),:,:) = 2;
% simroiimg((1+round(2 * dim(1)/3)):dim(1),:,:) = 3;
% 
% 
% simimg = zeros([dim length(b)]);
% 
%     
% for x=1:dim(1)
%     for y=1:dim(2)
%         for z=1:dim(3)
%             roi = simroiimg(x,y,z);
%             
%             for i=1:length(b)
%                 S(i) = f(roi) * exp(-b(i) * Dv(roi)) ...
%                     + (1 - f(roi)) * exp(-b(i) * D(roi));  
%                                 
%             end
%             
%             S = add_noise(S,SNR,'rician');
%                                          
%             simimg(x,y,z,:) = S;
%         end
%     end
% end
% 
% 
% % fit inspect segmentation model
% mask = ones(dim);
% 
% %inspect segmentation version
% 
% ILT_options = default_ILT_options('D');
% 
% inspect_seg_options.ILT = ILT_options;
% 
% %the number of clusters is known from the simulation
% inspect_seg_options.EM.n_clusters = nroi;
% inspect_seg_options.EM.n_steps = 10;
% inspect_seg_options.EM.init = 'kmeans';
% inspect_seg_options.save = 0;
% %inspect_options.save_path = [paperpath '/simulations/'];
% 
% 
% d_output = inspect_seg(simimg,gradechoinv,mask,'D',inspect_seg_options);
% 
% 
% 
% 
% 
% % %% T1-T2-D example
% % 
% % dim = [10 10 1];
% % 
% % gradechoinv = load('~/Dropbox/challenge/parameters.txt');
% % gradechoinv(:,7) = 7500;
% % %these are the other way around!
% % gradechoinvtemp = gradechoinv;
% % gradechoinv(:,5) = gradechoinvtemp(:,6);
% % gradechoinv(:,6) = gradechoinvtemp(:,5);
% % 
% % d = [0.002 0.001 0.003 0.002];
% % k = [0.7 1 0 0];
% % t2 = [100 100 150 150];
% % t1 = [500 1000 1500 2000];
% % 
% % SNR = 10000;
% % 
% % nroi = 4;
% % 
% % simroiimg = ones(dim);
% % 
% % simroiimg(1:2,:,:) = 1;
% % simroiimg(3:4,:,:) = 2;
% % simroiimg(5:7,:,:) = 3;
% % simroiimg(8:10,:,:) = 4;
% % 
% % 
% % simimg = zeros([dim size(gradechoinv,1)]);
% % 
% % for x=1:dim(1)
% %     for y=1:dim(2)
% %         for z=1:dim(3)
% %             roi = simroiimg(x,y,z);
% %             
% %             %params = [d(roi) t2(roi) t1(roi)];
% %             params = [d(roi) k(roi) t2(roi) t1(roi)];
% %             
% %             %S = KernelDT2T1(params,gradechoinv);
% %             S = KernelDKT2T1(params,gradechoinv);
% %     
% %             S = add_noise(S,SNR,'rician');
% %                                          
% %             simimg(x,y,z,:) = S;
% %         end
% %     end
% % end
% % 
% % 
% % % fit inspect segmentation model
% % mask = ones(dim);
% % 
% % clear inspect_seg_options
% % 
% % %the number of clusters is known from the simulation
% % inspect_seg_options.nclus = 4;
% % inspect_seg_options.nstep = 5;
% % 
% % inspect_seg_options.maxiter = 10;
% % 
% % 
% % dkt2t1_output = inspect_seg(simimg,gradechoinv,mask,'DKT2T1',inspect_seg_options);
% % 
% 
% 
% 
% 
% 
