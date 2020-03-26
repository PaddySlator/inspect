addpath(genpath('~/MERA-ext/'))

datadir = '~/challenge/';

img = [datadir 'cdmri11_r.nii'];
mask = [datadir 'cdmri11_mask.nii'];

imgs = {'cdmri12.nii','cdmri13.nii','cdmri14.nii','cdmri15.nii','cdmri11.nii'};
masks = {'cdmri12_mask.nii','cdmri13_mask.nii','cdmri14_mask.nii','cdmri15_mask.nii','cdmri11_mask.nii'};


gradechoinv = importdata([datadir 'parameters_new.txt']);
TR = importdata([datadir 'tr_MB2.txt']);

gradechoinvtemp = gradechoinv;
gradechoinv(:,5) = gradechoinvtemp(:,6);
gradechoinv(:,6) = gradechoinvtemp(:,5);
gradechoinv(:,7) = TR;


ILT_options = default_options('DT2T1');

inspect_options.ILT = ILT_options;
inspect_options.ILT.alpha = 0.001;

%the number of clusters is known from the simulation
inspect_options.EM.n_clusters = 6;
inspect_options.EM.n_steps = 5;
inspect_options.EM.init = 'random';
inspect_options.save = 1;

for i=1:length(imgs)
	img=[datadir imgs{i}];
	mask=[datadir masks{i}];
	dt2t1_output_brain = inspect_seg(img,gradechoinv,mask,inspect_options);
end

