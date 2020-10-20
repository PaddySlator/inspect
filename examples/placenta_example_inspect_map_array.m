%%% THE IMPORTANT ARRAY JOB STUFF
%get task id (i.e. index) of this job within the array job
job_id_string=getenv('SGE_TASK_ID')

if isempty(job_id_string)
   job_id=-1; 
else
   job_id=str2double(job_id_string);
end
 
%%add everything to path
addpath(genpath('dwi_placenta'))
addpath(genpath('downloaded_packages'))
addpath(genpath('MERA-ext'))


%% get the filenames of the data to fit to

%full path of the directory containing images, masks, protocol files, etc 
imgpath = '/Users/paddyslator/data/t2sdiff/';
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


%get the filenames of all the scans with include = 1
l=1;
scandir = cell(sum(include),1);
imgfilepath = cell(sum(include),1);
maskfilepath = cell(sum(include),1);
scanname = cell(sum(include),1);
gradechoinvfilepath = cell(sum(include),1);

for i=1:length(pipid)
    if include(i)
        scandir{l} = [imgpath pipid{i}];
        
        imginfo = dir([imgpath pipid{i} '/*alle.nii*']); 
        maskinfo = dir([imgpath pipid{i} '/*placenta_and_uterine_wall_mask.nii*']);
        gradechoinvinfo = dir([imgpath pipid{i} '/grad_echo_inv.txt']);
        
        scanname{l}  = pipid{i};
        imgfilepath{l} = [imgpath pipid{i} '/' imginfo.name];
        maskfilepath{l} = [imgpath pipid{i} '/' maskinfo.name];
        gradechoinvfilepath{l} = [imgpath pipid{i} '/' gradechoinvinfo.name];
        
        l=l+1;
    end
end










%% inspect options

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


%number of compartments (will be chosen later depending on the job id)
%inspect_options.n_comp = 2 ;

inspect_options.maxiter = 20 ;
inspect_options.init = 'random';
inspect_options.init = 'kmeans';
inspect_options.init = 'meanspectrum';
%inspect_options.init = 'user';

inspect_options.onhill = 0;
inspect_options.fmincon =  1;
inspect_options.updateF = 0;

inspect_options.parallel = 1;

inspect_options.hill.nsteps = 20;
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
inspect_options.dirname = 'placenta_fits';

inspect_options.save_path = '/Users/paddyslator/Documents/PlacentaDocs/papers/inspect_map/';
inspect_options.save_path = '/home/pslator/IPMI2019/';


%% run inspect map 

%batch job with different number of compartments and different scans

%max/min compartment numbers
maxcomp = 10;
mincomp = 2;


%get unique diffusion scan and model to fit 

n_data = length(imgfilepath);
n_totalcomps = maxcomp - mincomp;
n_jobs = n_data * n_totalcomps;

grid_combs=ndgrid(1:n_data,mincomp:maxcomp);
vec_combs(1,:)=grid_combs(:);
grid_combs=(ndgrid(mincomp:maxcomp, 1:n_data))';
vec_combs(2,:)=grid_combs(:);

%get unique indices (i.e. scan and number of compartments) from the job_id  
n_data
job_id
%scan to use
i=vec_combs(1,job_id)
%number of compartments 
j=vec_combs(2,job_id)

inspect_options.n_comp = j;
inspect_options.scan_names = scanname(i);

placentainspectmap = inspect_map(imgfilepath{i}, gradechoinvfilepath{i},maskfilepath{i},inspect_options);



        






