%%add everything to path
addpath(genpath('dwi_placenta'))
addpath(genpath('downloaded_packages'))
addpath(genpath('MERA-ext'))


%% load data

%full path of the directory containing images, masks, protocol files, etc 
imgpath = '/Users/paddyslator/Dropbox/placentaJhu/t2sdiff/';
imgpath = '/home/pslator/data/t2sdiff/';

paperpath = '/Users/paddyslator/Documents/PlacentaDocs/papers/IPMI2019/';
paperpath = '/home/pslator/IPMI2019/';


[~,~,raw] = xlsread([paperpath 'ipmipipdetails']);

%mask - fit only to voxels inside the mask
maskname = 'placenta_and_uterine_wall_mask';


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


%% put all the images together into a single structure
% fit to mutliple images - try to get two spectra representing "normal" and
%"bad" tissue
%concatenate the images together
if exist('t2diff','var')
    clear imgs masks gradechoinv
end

l=1;
for i=1:length(pipid)
    if include(i)
        imgs{l} = t2diff.(pipid{i}).normt2diffimg.img;
        masks{l} = t2diff.(pipid{i}).masks.(maskname).img;
        
        inspect_options.scan_names{l} = pipid{i};

        
        %they all have the same gradechoinv
        gradechoinv = t2diff.(pipid{i}).gradecho;
        
        l=l+1;               
    end
end




%% run inspect map 


%inspect options

%ILT options
ILT_options.Nk1=50;
ILT_options.Nk2=50;

ILT_options.Nk = [50 50];

ILT_options.mink1 = 2 * 10^-4;
ILT_options.maxk1 = 5000 * 10^-3;
ILT_options.mink2 = 10 * 10^-3;
ILT_options.maxk2 = 150 * 10^-3;

ILT_options.mink = [2*10^-4  10*10^-3];
ILT_options.maxk = [5000*10^-3  150*10^-3];


ILT_options.reg = 0;
%regularisation parameter
ILT_options.alpha = 0;

%regularisation for the mean fit
inspect_options.ILT_mean = ILT_options;
inspect_options.ILT_mean.alpha = 0.01;

%regularisation for the inspect map 
inspect_options.ILT = ILT_options;


%number of compartments 
inspect_options.n_comp = 3 ;

inspect_options.maxiter = 5 ;
inspect_options.init = 'random';
inspect_options.init = 'kmeans';
inspect_options.init = 'meanspectrum';
%inspect_options.init = 'user';

inspect_options.onhill = 0;
inspect_options.fmincon =  1;
inspect_options.updateF = 0;

inspect_options.parallel = 1;

inspect_options.hill.nsteps = 100;
inspect_options.hill.stepsize = 0.05;
inspect_options.hill.adaptwindow = 10;

inspect_options.hill.stepsizeadjustfactor = 5;

inspect_options.hill.adjustlowertol = 10^-2;
inspect_options.hill.adjustuppertol = 10;

inspect_options.weightstol = 10^-3;

inspect_options.onF = 1;
inspect_options.onweights = 1;

inspect_options.relabel = 1;

inspect_options.save = 1;
inspect_options.dirname = 'placenta_fits_allimgs';
inspect_options.dirname = 'placenta_fits';
inspect_options.save_path = '/Users/paddyslator/Documents/PlacentaDocs/papers/inspect_map/';
inspect_options.save_path = '/home/pslator/IPMI2019/'



%fit to single images 
%choose the image to fit to
%scanname = 'pip0111';

%loop over all images 
%choose these images

pipid = {'pip012002','pip0120','pip0108'};
ncomp = [3 4 5];
include = [1 1 1];

for i=1:length(pipid)
    if include(i)
            
        %load the placenta gradechoinv file
        %gradechoinv = load('/Users/paddyslator/Dropbox/placentaJhu/t2sdiff/pip0111/grad_echo_inv.txt');
        
        gradechoinv = load([imgpath (pipid{i}) '/grad_echo_inv.txt']);

        nmeas = size(gradechoinv,1);
        
        img = t2diff.( pipid{i} ).normt2diffimg.img;
        mask = t2diff.( pipid{i} ).masks.(maskname).img;
        
        %loop over a number of different compartments
        %maxcomp = 2;
        %mincomp = 10;
        %for j=mincomp:maxcomp
            %inspect_options.n_comp = j;
	    inspect_options.n_comp = ncomp(i)
	    inspect_options.scan_names = pipid{i}
            placentainspectmap = inspect_map(img, gradechoinv,mask,inspect_options);
        %end
        
    end
end


% fit to all images at once
%placentainspectmap = inspect_map(imgs, gradechoinv,masks,inspect_options);









