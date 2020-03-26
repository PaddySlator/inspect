%% load data

%full path of the directory containing images, masks, protocol files, etc 
imgpath = '/Users/paddyslator/Dropbox/placentaJhu/t2sdiff/';

paperpath = '/Users/paddyslator/Documents/PlacentaDocs/papers/IPMI2019/';

[~,~,raw] = xlsread([paperpath 'ipmipipdetails']);


%get the scans to include
includecell = strfind(raw(1,:),'include');
includecol = find(not(cellfun('isempty',includecell)));
include = cell2mat(raw(2:end,includecol));

%get the pip ids
pipcell = strfind(raw(1,:),'pipid');
pipcol = find(not(cellfun('isempty',pipcell)));
pipid = raw(2:end,pipcol);

%get the cohort
cohortcell = strfind(raw(1,:),'cohort');
cohortcol = find(not(cellfun('isempty',cohortcell)));
cohort = raw(2:end,cohortcol);


loadon = 1;
if loadon 
    %load all the scans with include = 1
    for i=1:length(pipid)
        if include(i)
            t2diff.(pipid{i}) = load_T2MEdiff([imgpath pipid{i}]);
            disp(['loaded scan ' pipid{i}])
        end        
    end
end



%% run inspect

%mask - fit only to voxels inside the mask
maskname = 'placenta_and_uterine_wall_mask';

% %just fit to the 10th slice for testing
% for i=1:length(pipid)
%     if include(i)
%         newmask = zeros(size(t2diff.(pipid{i}).masks.(maskname).img));
%         newmask(:,:,10) = t2diff.(pipid{i}).masks.(maskname).img(:,:,10);
%         t2diff.(pipid{i}).masks.(maskname).img = newmask;
%     end
% end



%inspect options

ILT_options.Nk1=50;
ILT_options.Nk2=50;

ILT_options.mink1 = 2 * 10^-4;
ILT_options.maxk1 = 5000 * 10^-3;
ILT_options.mink2 = 10 * 10^-3;
ILT_options.maxk2 = 150 * 10^-3;
ILT_options.alpha = 0.01;


inspect_options.ILT = ILT_options;
inspect_options.EM.n_clusters = 3;
inspect_options.EM.n_steps = 50;
inspect_options.EM.init = 'random';


inspect_options.save=1;

inspect_options.save_path = paperpath;




% fit to mutliple images - try to get two spectra representing "normal" and
%"bad" tissue
%concatenate the images together
if exist('t2diff','var')
    clear img mask gradechoinv
end

l=1;
for i=1:length(pipid)
    if include(i)
        img{l} = t2diff.(pipid{i}).normt2diffimg.img;
        mask{l} = t2diff.(pipid{i}).masks.(maskname).img;
        
        inspect_options.scan_names{l} = pipid{i};

        
        %they all have the same gradechoinv
        gradechoinv = t2diff.(pipid{i}).gradecho;
        
        l=l+1;               
    end
end

%hopefully this will save some memory!
clear t2diff



%% fit to all images individually

clusters = 2:10;

l=1;
for i=1:length(pipid)
    inspect_options.save_path = [paperpath pipid{i}];
    if include(i)
        for j=1:length(clusters)
            
            inspect_options.EM.n_clusters = clusters(j);
            
            [~,inspectoutputsummary.(pipid{i}){j}] = inspect_seg(img{l},...
                gradechoinv,...
                mask{l},...
                inspect_options);               
        end
        l=l+1;
        disp(['done pip ' pipid{i}])
    end
end


%% get some stuff out from the fits
l=1;
for i=1:length(pipid)
    if include(i)
        for j=1:length(clusters)
            BIC(l,j) = inspectoutputsummary.(pipid{i}){j}.BIC;
            AIC(l,j) = inspectoutputsummary.(pipid{i}){j}.AIC;
        end
        l=l+1;
    end
end
           


%% voxelwise fit




%% run on all images at once


clusters = 2:10;

inspect_options.save=1;

l=1;
for j=clusters(1):clusters(end)
    inspect_options.EM.n_clusters = j;
    inspect_options.save_path = [paperpath 'allimg_' num2str(j) '_clusters'];
    [~,inspectoutputmultisummary.(['clusters' num2str(j)])] = inspect_seg(img,gradechoinv,mask,inspect_options);
    
end



%% get BIC out
for j=clusters
    BICalldata(j) = inspectoutputmultisummary.(['clusters' num2str(j)]).BIC;
    AICalldata(j) = inspectoutputmultisummary.(['clusters' num2str(j)]).AIC;
end





return

%% fit to a single image
img = double(t2diff.t2diffimg.img);
mask = double(t2diff.masks.(maskname).img);


newmask = zeros(size(mask));
newmask(:,:,10) = mask(:,:,10);
mask = newmask;

gradechoinv = double(t2diff.gradecho);




%give the correct answer as a start point
%inspect_options.EM.initweights = weightsimg;

%fit to each image
inspectoutput = cell(size(imgpath));
for i=1:length(imgpath)
    inspectoutput{i} = inspects(t2diff{i}.t2diffimg.img,...
        t2diff{i}.gradecho,...
        t2diff{i}.masks.(maskname).img,...
        inspect_options);
end




%% load MERA voxelwise fits


imgpath = '/Users/paddyslator/Dropbox/placentaJhu/t2sdiff/';

if loadon 
    %load all the scans with include = 1
    for i=1:length(pipid)
        if include(i)
            try
                MERAvoxel.(pipid{i}) = load([imgpath pipid{i} '/placenta_and_uterine_wall_mask_single_slice_MERA_fit.mat']);
                
                maskforMERA.(pipid{i}) = load_untouch_nii([imgpath pipid{i} '/placenta_and_uterine_wall_mask.nii.gz']);
                disp(['loaded scan ' pipid{i}])
            catch
                disp(['couldnt load scan ' pipid{i}])
            end
        end        
    end
end


%% 

w1bounds = [0 25; 0 25; 0 25; 25 200 ;25 200 ;25 200 ;200 1000; 200 1000; 200 1000];
w2bounds = [0 0.03; 0.03 0.07;0.07 1;0 0.03; 0.03 0.07;0.07 1;0 0.03; 0.03 0.07;0.07 1];

w1bounds= [0 20; 20 1000; 20 1000;0 20];
w2bounds= [0 0.075; 0.075 0.2;0 0.075; 0.075 0.2];

w1bounds= [0 1000; 0 1000];
w2bounds= [0 0.075; 0.075 0.2];
    
w1bounds= [0 10; 0 10; 0 10; 10 1000; 10 1000; 10 1000];
w2bounds= [0 0.025; 0.025 0.075; 0.075 1; 0 0.025; 0.025 0.075; 0.075 1];

%w1bounds = [0 10;10 1000;0 10;10 1000;0 10; 10 1000];
%w2bounds = [0 0.025;0 0.025;0.025 0.075;0.025 0.075;0.075 0.2;0.075 0.2];


for i=1:length(pipid)
    if include(i)
        %use the mask as a template
        Vf = zeros([size(maskforMERA.(pipid{i}).img) size(w1bounds,1)]);       
        for j=1:(length(MERAvoxel.(pipid{i}).MERA_fit)-1)
            voxind = MERAvoxel.(pipid{i}).MERA_fit{j}.voxel;
            Vf(voxind(1),voxind(2),voxind(3),:) = integrate_spectrum_2D(MERAvoxel.(pipid{i}).MERA_fit{j},w1bounds,w2bounds,0);
        end
        
        Vfnii = make_nii(Vf);
        save_nii(Vfnii,[imgpath pipid{i} '/placenta_and_uterine_wall_mask_single_slice_MERA_fv_map_' 'adc_0_10_1000_t2_0_025_075_1.nii.gz'])
    end
end









