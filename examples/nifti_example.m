img_filename = 'example_dataset/pip0130_trunc.nii.gz';
gradechoinv_filename = 'example_dataset/grad_echo_inv.txt';
mask_filename = 'example_dataset/placenta_mask_trunc.nii.gz';

options.average_fit = 1;
%warning: if you set this to 1 it will take ages (probably hours)!!! 
options.voxelwise_fit = 0;


% choose MERA fitting options
% fitting.twoD = 'y';
% analysis.interactive = 'n'
% fitting.regtyp = 'me';
% %DIFFUSION
% fitting.rangeT = [2,5000] * 10^-3;
% %RELAXOMETRY
% fitting.rangeT2 = [10,150] * 10^-3 ;
% fitting.numbergauss = 3;
% fitting.regadj = 'manual';
% analysis.graph = 'y';
% fitting.regweight = 0.002;
% fitting.numberT = 50;
% fitting.numberT2 = 50;
% analysis.extract = 'auto'
% analysis.graph = 'n';
% disp('using default MERA fitting options')
% MERA_options.fitting = fitting;
% MERA_options.analysis = analysis;
% MERA_fit = MERA_fit_nifti(img_filename,GradEchoInv_filename,mask_filename,options,MERA_options);

    
%or ... just give 4 arguments - so will use the default MERA fitting options
MERA_fit = MERA_fit_nifti(img_filename,gradechoinv_filename,mask_filename,options);

