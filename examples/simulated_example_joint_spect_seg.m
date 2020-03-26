
%simulate a simple image where each voxel has an associated spectrum 
dim = [10 10 10];

%load a placenta gradechoinv file
gradechoinv = load('/Users/paddyslator/Dropbox/t2sdiff/pip0111/grad_echo_inv.txt');
nmeas = size(gradechoinv,1);

%simroiimg = zeros(dim);

%number of ROIs - each with a different associated spectrum
nroi = 3;

%
%simweightsimg = eps + zeros([dim 3]);


% %associated spectrum values for each ROI
% %[T2 D f]
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




%% do a load of single compartment spectra!

SNR=10^6;

simoptions.SNR = SNR;
simoptions.noisetype = 'rician';

disp(['SNR SET TO ' num2str(SNR)])

%experiments with ADC fixed 
T2base = [0.02 0.05 0.08];
ADCbase = [0.002 0.02 0.2];

adcscalings = [1 1 1 1];
t2scalings = [1.25 1.5 1.75 2];



ntest = length(adcscalings)*length(ADCbase)*length(T2base);

%spectparams=cell(ntest,1);

%simimg = cell(1);
%simweightsimg = cell(1);

for i=1:length(T2base)
    for j=1:length(ADCbase)
        for k=1:length(t2scalings)
            %first ROI/spectra - always the same
            spectparams.t2{i}{j}{k}{1}.T2 = T2base(i);
            spectparams.t2{i}{j}{k}{1}.D = ADCbase(j);
            spectparams.t2{i}{j}{k}{1}.f = 1;
            
            %second ROI/spectra - moves away in T2* dim
            spectparams.t2{i}{j}{k}{2}.T2 = t2scalings(k)*T2base(i);
            spectparams.t2{i}{j}{k}{2}.D = ADCbase(j);
            spectparams.t2{i}{j}{k}{2}.f = 1;
            
            %third ROI/spectra - moves away in ADC dim
            %spectparams{i}{3}.T2 = T2base;
            %spectparams{i}{3}.D =  adcscalings(i)*ADCbase;
            %spectparams{i}{3}.f = 1;
            
            %fourth ROI/spectra - moves away in T2* and ADC dim
            %spectparams{i}{4}.T2 = t2scalings(i)*T2base;
            %spectparams{i}{4}.D =  adcscalings(i)*ADCbase;
            %spectparams{i}{4}.f = 1;
        end
    end
end

%make a checkerboard thing
simroiimg = zeros(dim);

simroiimg(1:5,:,:) = 1;
simroiimg(6:10,:,:) = 2;
%simroiimg(6:10,:,:) = 3;
%simroiimg(6:10,6:10,:) = 4;


for i=1:length(T2base)
    for j=1:length(ADCbase)
        for k=1:length(t2scalings)
            %now simulate the experiment
            for x=1:dim(1)
                for y=1:dim(2)
                    for z=1:dim(3)
                        %sample the ROI
                        roi = simroiimg(x,y,z);
                                                                    
                        S = simulate_multiexp_signal(spectparams.t2{i}{j}{k}{roi},gradechoinv,simoptions);                        
                        
                        simimg.t2{i}{j}{k}(x,y,z,:) = S;
                        simweightsimg.t2{i}{j}{k}(x,y,z,roi) = 1 - eps;
                    end
                end
            end
        end
    end
end



%SIMULATE CHANGING THE ADC

adcscalings = [1.25 1.5 1.75 2];
t2scalings = [1 1 1 1];

for i=1:length(T2base)
    for j=1:length(ADCbase)
        for k=1:length(adcscalings)
            %first ROI/spectra - always the same
            spectparams.adc{i}{j}{k}{1}.T2 = T2base(i);
            spectparams.adc{i}{j}{k}{1}.D = ADCbase(j);
            spectparams.adc{i}{j}{k}{1}.f = 1;
            
            %second ROI/spectra - moves away in ADC dim
            spectparams.adc{i}{j}{k}{2}.T2 = T2base(i);
            spectparams.adc{i}{j}{k}{2}.D = adcscalings(k)*ADCbase(j);
            spectparams.adc{i}{j}{k}{2}.f = 1;                        
        end
    end
end


for i=1:length(T2base)
    for j=1:length(ADCbase)
        for k=1:length(adcscalings)
            %now simulate the experiment
            for x=1:dim(1)
                for y=1:dim(2)
                    for z=1:dim(3)
                        %sample the ROI
                        roi = simroiimg(x,y,z);
                        
                       
                        
                        S = simulate_multiexp_signal(spectparams.adc{i}{j}{k}{roi},gradechoinv,SNR);
                        
                        
                        simimg.adc{i}{j}{k}(x,y,z,:) = S;
                        simweightsimg.adc{i}{j}{k}(x,y,z,roi) = 1 - eps;
                    end
                end
            end
        end
    end
end




%% fit the spectrum to each voxel individually 

output=cell(ntest,1);

ILT_options.Nk1=50;
ILT_options.Nk2=50;

ILT_options.mink1 = 2 * 10^-4;
ILT_options.maxk1 = 5000 * 10^-3;
ILT_options.mink2 = 10 * 10^-3;
ILT_options.maxk2 = 150 * 10^-3;

ILT_options.alpha = 0.01;

%i=2;j=2;k=1;
    ... for i=1:length(T2base)
        ...   for j=1:length(ADCbase)
        ...    for k=1:length(adcscalings)

k=4;
for change={'adc','t2'}
    for i = 1:3
        for j = 1:3
            %output{i}{j}{k} = cell(dim);
            for x=1:dim(1)
                for y=1:dim(2)
                    for z=1 %only do one slice! only need this for the images                       
                        sig = squeeze(simimg.(change{1}){i}{j}{k}(x,y,z,:));
                        %output.(change{1}){i}{j}{k}{x,y,z} = ILT_2D(sig, gradechoinv,ILT_options);
                        
                        %don't store this because it becomes massive
                        output = ILT_2D(sig, gradechoinv,ILT_options);
                        %calculate vf map
                        if strcmp(change{1},'adc')
                            midpoint = (spectparams.(change{1}){i}{j}{k}{1}.D + spectparams.(change{1}){i}{j}{k}{2}.D)/2;
                            boundadc = [0 midpoint;midpoint 1000];
                            boundt2 = [0 1;0 1];
                        else
                            midpoint = (spectparams.(change{1}){i}{j}{k}{1}.T2 + spectparams.(change{1}){i}{j}{k}{2}.T2)/2;
                            boundt2 = [0 midpoint;midpoint 1];
                            boundadc = [0 1000;0 1000];
                        end
                        
                        Vf.(change{1}){i}{j}{k}(x,y,z,:) = integrate_spectrum_2D(output,boundadc,boundt2);                        
                                           
                    end
                end
            end
            disp(change{1})
            disp(['i = ' num2str(i)])
            disp(['j = ' num2str(j)])
        end
    end
end

%correctvoxvoxelwise{i}{j}{k} = sum(sum(sum()));





%% test calculate mean signal a bit




%% do inspects!
mask = ones(dim);

inspect_options.ILT = ILT_options;
inspect_options.EM.n_clusters = 2;
inspect_options.EM.n_steps = 10;
inspect_options.EM.init = 'random';

paperpath = '/Users/paddyslator/Documents/PlacentaDocs/papers/IPMI2019/';

saveon=0
inspect_options.save = saveon;
inspect_options.save_path = [paperpath '/simulations/'];

inspect_options.scan_names{1} = ['SNR_' num2str(SNR) '_' num2str(prod(dim)) '_voxels' ];

%give the correct answer as a start point
%inspect_options.EM.initweights = weightsimg;

for change={'adc','t2'}
    siminspectoutput.(change{1})=cell(ntest,1);
    for i=1:length(T2base)
        for j=1:length(ADCbase)
            for k=1:length(t2scalings)
                siminspectoutput.(change{1}){i}{j}{k} = inspect_seg(simimg.(change{1}){i}{j}{k},gradechoinv,mask,inspect_options);
            end
        end
    end
end

%% check if they need swapping

indices=[1:5;6:10];

for change={'adc','t2'}
    for i=1:length(T2base)
        for j=1:length(ADCbase)
            for k=1:length(t2scalings)
                meanroi(1) = mean(mean(mean(siminspectoutput.(change{1}){i}{j}{k}.MLroi{1}(1:5,:,:))));
                meanroi(2) = mean(mean(mean(siminspectoutput.(change{1}){i}{j}{k}.MLroi{1}(6:10,:,:))));
                
                if meanroi(2)<meanroi(1)
                    tempMLroi = siminspectoutput.(change{1}){i}{j}{k}.MLroi{1};
                    siminspectoutput.(change{1}){i}{j}{k}.MLroi{1}(indices(2,:),:,:) = tempMLroi(indices(1,:),:,:);
                    siminspectoutput.(change{1}){i}{j}{k}.MLroi{1}(indices(1,:),:,:) = tempMLroi(indices(2,:),:,:);
                end
            end
        end
    end
end


%% get some stuff out and save
for change={'adc','t2'}
    for i=1:length(T2base)
        for j=1:length(ADCbase)
            for k=1:length(t2scalings)
                MLroi.(change{1}){i}{j}{k} =  siminspectoutput.(change{1}){i}{j}{k}.MLroi{1};
                imgweights.(change{1}){i}{j}{k} =  siminspectoutput.(change{1}){i}{j}{k}.imgweights{1};
                BIC.(change{1}){i}{j}{k} =  siminspectoutput.(change{1}){i}{j}{k}.AIC;
                AIC.(change{1}){i}{j}{k} =  siminspectoutput.(change{1}){i}{j}{k}.BIC;
                sim_ILT_output.(change{1}){i}{j}{k} =  siminspectoutput.(change{1}){i}{j}{k}.ILT_output;
            end
        end
    end
end
if saveon
    save('INSPECT_SIM.mat','MLroi','imgweights','BIC','AIC')
end

%% plot loads of voxelwise vs thingy comparisons!
figure;
k=4;
l=1;
ncol=1;
nrow=9;
for change={'adc','t2'}
    figure;
    l=1;
    for i=1:length(T2base)
        for j=1:length(ADCbase)            
            subtightplot(ncol,nrow,l);
            imagesc(Vf.(change{1}){i}{j}{k}(:,:,1,1))
            axis off
            axis square
            l=l+1;
            
            title(['(' num2str(spectparams.t2{i}{j}{k}{1}.T2) ',' num2str(spectparams.t2{i}{j}{k}{1}.D) '); ' '(' num2str(spectparams.t2{i}{j}{k}{2}.T2) ',' num2str(spectparams.t2{i}{j}{k}{2}.D) ')']) 
            if l==2
              % title('Voxelwise Map') 
            end
        end
    end
    set(gcf,'Position',[56        1144        1256         511])
        
    print_to_formats([figdir '/ROI_VOXELWISE_change' change{1}],{'-depsc','-dpng','fig'})

    figure;
    l=1;
    for i=1:length(T2base)
        for j=1:length(ADCbase)            
            subtightplot(ncol,nrow,l);
            imagesc(imgweights.(change{1}){i}{j}{k}(:,:,1))
            axis off
            axis square
            l=l+1;
            
            title(['(' num2str(spectparams.t2{i}{j}{k}{1}.T2) ',' num2str(spectparams.t2{i}{j}{k}{1}.D) '); ' '(' num2str(spectparams.t2{i}{j}{k}{2}.T2) ',' num2str(spectparams.t2{i}{j}{k}{2}.D) ')']) 
            if l==2
               % title('InSpect Map')
            end
        end
    end
    set(gcf,'Position',[56        1144        1256         511])

    print_to_formats([figdir '/ROI_INSPECT_change' change{1}],{'-depsc','-dpng','fig'})

end




%% get the proportions of voxels correctly classified for each experiment
for change={'adc','t2'}
    for i=1:length(T2base)
        for j=1:length(ADCbase)
            for k=1:length(t2scalings)
                correctvox.(change{1})(i,j,k) = sum(sum(sum(siminspectoutput.(change{1}){i}{j}{k}.MLroi{1} == simroiimg)));
            end
        end
    end
    
    if strcmp(change{1},'adc')
        for n=1:length(adcscalings)
            adcscalingcellstring{n} = num2str(adcscalings(n));
        end
    else
        for n=1:length(t2scalings)
            t2scalingcellstring{n} = num2str(t2scalings(n));
        end
    end
end

%% plot these proportions
gap = 0.03;
marg_h = 0.12;
marg_v = 0.05;

for change={'adc','t2'}
    figure;
    for i=1:3
        subtightplot(1,3,i,gap,marg_h,marg_v);hold on;
        %SAME T2*/ADC BASE IN EACH SUBPLOT T2* is first dim, ADC is second
        if strcmp(change{1},'adc')
            bar(squeeze(correctvox.(change{1})(:,i,:))'./prod(dim))
        else
            bar(squeeze(correctvox.(change{1})(i,:,:))'./prod(dim))
        end
        if strcmp(change{1},'adc')
            title(['First component''s ADC = ' num2str(ADCbase(i)) ' mm^2/s'],'FontSize',18)
        else
            title(['First component''s T2^* = ' num2str(T2base(i)) ' s'],'FontSize',18)
        end
        
        if strcmp(change{1},'adc')
            set(gca,'XTick',(1:length(adcscalings)))
            set(gca,'XTickLabel',adcscalingcellstring)
            xlabel('ADC Ratio between spectral components','FontSize',14)
        else
            set(gca,'XTick',(1:length(t2scalings)))
            set(gca,'XTickLabel',adcscalingcellstring)
            xlabel('T2^* Ratio between spectral components','FontSize',14)
        end
        ylim([0 1])
        if i>1
            set(gca,'yticklabel',{[]})
        else
            ylabel('Proportion of voxels clustered correctly','FontSize',14)
        end
    end
    
    if strcmp(change{1},'adc')
        for j=1:length(T2base)
            legstring{j} = ['T2^* = ' num2str(T2base(j)) ' mm^2/s'];
        end
    else
        for j=1:length(ADCbase)
            legstring{j} = ['ADC = ' num2str(ADCbase(j)) ' mm^2/s'];
        end
    end
    
    legend(legstring,'FontSize',10,'Location','southwest')
    set(gcf,'Position',[ 288        1348        1161         340])
    
    print_to_formats([figdir '/SIM_VARY_' change{1}],{'-depsc','-dpng','fig'})
    
end

%plot the positions of the two peaks




%% plot the roi
figdir = '/Users/paddyslator/Documents/PlacentaDocs/papers/IPMI2019/figures';

nroi = inspect_options.EM.n_clusters;


i=1;j=3;k=1;

ij=[3 2 1;1 2 3];

gap = [0.03 0.02];
marg_height = [0.07 0.05];
marg_width = [0.07 0.05];

for change={'adc','t2'}
    figure;
    l=1;
    plotdim=[3 2];
    for ij=[3 2 1;1 2 3]
        %plot the spectra
        subtightplot(plotdim(1),plotdim(2),l,gap,marg_width,marg_height);hold on;l=l+1;
        
        %need to match up the inferred clusters with the simulated clusters - could
        %do this automatically quite easily
        %roiordering = [1 2 3 4];
        
        
        for m=1:nroi
            %subplot(1,nroi,m);hold on;
            
            contour(siminspectoutput.(change{1}){ij(1)}{ij(2)}{k}.ILT_output{m}.w2,...
                siminspectoutput.(change{1}){ij(1)}{ij(2)}{k}.ILT_output{m}.w1,...
                siminspectoutput.(change{1}){ij(1)}{ij(2)}{k}.ILT_output{m}.F)
            set(gca,'xscale','log')
            set(gca,'yscale','log')
            if l==6
                xlabel('T2^* (s)','FontSize',16)
            end
            if l==4
                ylabel('ADC (10^{-3} mm^2 s^{-1})','FontSize',16)
            end
            
            if l==2
                title(['InSpect Spectra'],'FontSize',20);
            end
            
            colors={'r','k'};
            %plot simulation peaks
            %for n=1:length(spectparams{i}{j}{k}{m})
            plot(spectparams.(change{1}){ij(1)}{ij(1)}{k}{m}.T2,spectparams.(change{1}){ij(1)}{ij(2)}{k}{m}.D,[colors{m} 'x']);
            %end
        end
        
        subtightplot(plotdim(1),plotdim(2),l,gap,marg_width,marg_height);l=l+1;
        imagesc(siminspectoutput.(change{1}){ij(1)}{ij(2)}{k}.MLroi{1}(:,:,1))
        set(gca,'YDir','normal')
        axis off
        if l==3
            title('InSpect Maps','FontSize',20);
        end
        
        %print_to_formats([figdir '/inspectsimrois'],{'-depsc','fig'})
        set(gcf,'Position',[606   985   503   706])
    end
end



print_to_formats([figdir '/SIM_VARY_' change{1} '_image'],{'-depsc','-dpng','fig'})




%%
figure; hold on;
l=1;
for x=1:dim(1)
    for y=1:dim(2)
        subplot(dim(1),dim(2),sub2ind(dim,x,y));hold on;       
        contour(output{x,y}.w2,output{x,y}.w1*10^3,output{x,y}.F)
        set(gca,'xscale','log')
        set(gca,'yscale','log')      
        %plot simulation peaks
        for i=1:length(spectparams{simroiimg(x,y)}.T2)
            plot(spectparams{simroiimg(x,y)}.T2(i),spectparams{simroiimg(x,y)}.D(i)*10^3,'rx');
        end
        l=l+1;
    end
end
set (gcf,'Position',[216   881   921   689])
print_to_formats([figdir '/voxelwisesimpeaks'],{'-depsc','fig'})


%plot the ROIs
figure; hold on;
imagesc(simroiimg(:,:,1))
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
