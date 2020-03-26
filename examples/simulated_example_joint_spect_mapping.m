
%simulate a simple image where each voxel has multiple associated spectra - 
%and the signal in each voxel is a weighting of these spectra
dim = [4 4 4];


%load a placenta gradechoinv file
gradechoinv = load('/Users/paddyslator/Dropbox/placentaJhu/t2sdiff/pip0111/grad_echo_inv.txt');
nmeas = size(gradechoinv,1);

img = zeros([dim nmeas]);

%number of spectra - each voxel is a weight sum of some of these
nspect = 3;


%image of cluster/spectra memberships 
roiimg = randi(0:1,[dim nspect]);

%don't want any zeros
for x=1:dim(1)
    for y=1:dim(2)
        for z=1:dim(3)
            if sum(roiimg(x,y,z,:)) == 0
                roiimg(x,y,z,randi(nspect)) = 1;
            end
        end
    end
end

%image of weighting of each cluster
weightsimg = rand([dim nspect]);
%multiply by cluster memberships - make weight 0 if this cluster is not
%involved
weightsimg = weightsimg.*roiimg;
%normalise


for x=1:dim(1)
    for y=1:dim(2)
        for z=1:dim(3)
            weightsimg(x,y,z,:) = weightsimg(x,y,z,:)./sum(weightsimg(x,y,z,:));
        end
    end
end






%associated spectrum values for each cluster
%[T2 D f]
spectparams = cell(nspect,1);



spectparams{1}.T2 = [0.03];
spectparams{1}.D = [0.002];
spectparams{1}.f = [1];

spectparams{2}.T2 = [0.04];
spectparams{2}.D = [0.03];
spectparams{2}.f = [1];

spectparams{3}.T2 = [0.08];
spectparams{3}.D = [0.3];
spectparams{3}.f = [1];

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




%now simulate the experiment
for x=1:dim(1)
    for y=1:dim(2)
        for z=1:dim(3)
            %get the spectra memberships
            roi = roiimg(x,y,z,:);
            %get the spectra weights
            weights = weightsimg(x,y,z,:);
            
            SNR=1000;
            
            sigspect = zeros(nspect,size(gradechoinv,1));
            
            %signal component from each spectra
            for i=1:nspect
                sigspect(i,:) = weights(i)*simulate_multiexp_signal(spectparams{i},gradechoinv,SNR);
            end
            S=sum(sigspect);
            
            %normalise 
            b = gradechoinv(:,4);
            te = gradechoinv(:,5);
            
            S = S./(mean(S(b == 0 & te == min(te))));
            
            img(x,y,z,:) = S;
        end
    end
end





%% fit the spectrum to each voxel individually 
output = cell(dim);


ILT_options.Nk1=50;
ILT_options.Nk2=50;

ILT_options.mink1 = 2 * 10^-4;
ILT_options.maxk1 = 5000 * 10^-3;
ILT_options.mink2 = 10 * 10^-3;
ILT_options.maxk2 = 150 * 10^-3;

ILT_options.alpha = 0.01;
for x=1:dim(1)
    for y=1:dim(2)
        for z=1:dim(3)
            sig = squeeze(img(x,y,z,:));
            output{x,y,z} = ILT_2D(sig, gradechoinv,ILT_options);
        end
    end
end




%% do inspects!
mask = ones(dim);

inspect_options.ILT = ILT_options;
inspect_options.EM.n_clusters = 3;
inspect_options.EM.n_steps = 10;
inspect_options.EM.init = 'random';

paperpath = '/Users/paddyslator/Documents/PlacentaDocs/papers/IPMI2019/';

inspect_options.save = 0;
inspect_options.save_path = [paperpath '/simulations/'];

inspect_options.scan_names{1} = ['SNR_' num2str(SNR) '_' num2str(prod(dim)) '_voxels' ];

%give the correct answer as a start point
%inspect_options.EM.initweights = weightsimg;


inspectoutput = inspects_multi_img(img,gradechoinv,mask,inspect_options);


%% plot the roi
figdir = '/Users/paddyslator/Documents/PlacentaDocs/papers/IPMI2019/figures';


figure;
imagesc(inspectoutput.MLroi{1}(:,:,1))
set(gca,'YDir','normal')
axis off
title('INSPECTS Clustering')
print_to_formats([figdir '/inspectsimrois'],{'-depsc','fig'})


%plot the spectra
figure; hold on;

%need to match up the inferred clusters with the simulated clusters - could
%do this automatically quite easily
roiordering = [1 2 3];
for i=1:nspect
        subplot(1,nspect,i);hold on;       
        
        contour(inspectoutput.ILT_output{i}.w2,...
            inspectoutput.ILT_output{i}.w1,...
            inspectoutput.ILT_output{i}.F)
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        xlabel('T2^* (s)')
        ylabel('ADC (10^{-3} mm^2 s^{-1})')
        
        title(['Cluster ' num2str(i)]);
        
        %plot simulation peaks
        for j=1:length(spectparams{roiordering(i)}.T2)
            plot(spectparams{roiordering(i)}.T2(j),spectparams{roiordering(i)}.D(j),'rx');
        end    
end
set(gcf,'Position',[ 146         998        1054         278])

print_to_formats([figdir '/inspectsimpeaks'],{'-depsc','fig'})




%% plot the voxelwise fits
figure; hold on;
l=1;
for x=1:dim(1)
    for y=1:dim(2)
        subplot(dim(1),dim(2),sub2ind(dim,x,y));hold on;       
        contour(output{x,y}.w2,output{x,y}.w1*10^3,output{x,y}.F)
        set(gca,'xscale','log')
        set(gca,'yscale','log')      
        %plot simulation peaks
        %for i=1:length(spectparams{roiimg(x,y)}.T2)
        %    plot(spectparams{roiimg(x,y)}.T2(i),spectparams{roiimg(x,y)}.D(i)*10^3,'rx');
        %end
        l=l+1;
    end
end
set (gcf,'Position',[216   881   921   689])
print_to_formats([figdir '/voxelwisesimpeaks'],{'-depsc','fig'})


%plot the ROIs
figure; hold on;
imagesc(roiimg(:,:,1))
set(gca,'YDir','normal')
axis off
print_to_formats([figdir '/simrois'],{'-depsc','fig'})



%% save some stuff for plotting in python
% voxelwise spectra fits and simulated spectrum values
fitdir = '/Users/paddyslator/Documents/PlacentaDocs/papers/IPMI2019/simulations';

%only want the spectrum outputs
Fsim = cell(size(output));
w1sim = cell(size(output));
w2sim = cell(size(output));

for x=1:dim(1)
    for y=1:dim(2)
        for z=1:dim(3)
            Fsim{x,y,z} = output{x,y,z}.F;
            w1sim{x,y,z} = output{x,y,z}.w1;
            w2sim{x,y,z} = output{x,y,z}.w2;
        end
    end
end

saveon = 0;
if saveon
    save([fitdir '/voxelwise_fit.mat'],'spectparams','roiimg','img','Fsim','w1sim','w2sim')

    % inspects output
    save([fitdir '/inspects_fit.mat'],'spectparams','roiimg','img','inspectoutput')
end





%%
% 
% %fit with MERA 
% for x=1:dim(1)
%     for y=1:dim(2)
%         sig = squeeze(img(x,y,:));
%         MERA_data = make_MERA_data(sig,gradechoinv);
%         [MERA_out2D, MERA_fitting_out] = MERA_T2MEdiff_fit(MERA_data);
%         MERA_output{x,y} = MERA_out2D;    
%     end
% end
% 
% %%
% figure;hold on;
% contour(MERA_output{1,1}.T2,MERA_output{1,1}.T,MERA_output{1,1}.S)
% set(gca,'xscale','log')
% set(gca,'yscale','log')
