function output = inspect_map(img,gradechoinv,mask,kernel,options)

%inspect continuous mapping version
%assumes that the signal from each voxel is a weighted sum of 
%contributions from a set of canonicial basis spectra

% REQUIRED INPUTS:
%
% img - either a single image, a cell of multiple images, or file path to a
% single nifti image. Each image has dimension [dimx dimy dimz encodings]. 
% If using multiple images these must all have the same MR encodings (i.e. gradechoinv file). 
% For multiple images the model is fit to all images at once.
%
% gradechoinv - array, or path to a .txt file, of the MR acquisition 
% parameters.  Format: [gx gy gz b TE TI TR] where gx gy gz are the gradient 
% directions, b is the b-value, TE is the echo time, TI is the inversion 
% time, and TR is the repetition time. 
%
% mask - either a single binary image a cell of multiple binary images,
% or file path to a single nifti binary mask.
% Will be used to choose the voxels in which to fit the model. The mask
% must have dimension [dimx dimy dimz] - i.e. matching the first three
% dimensions of the corresponding image.
%
% kernel - string specifying the choice of kernel, which relates
% the MR acquisition parameters to the MR signal. For example:
% - 'DT2' for diffusion-T2, for which the kernel equation 
%   is exp(-b*D)*exp(-TE/t2). 
% - 'T1' for T1 inversion recovery, for which the kernel equation 
%   is abs(1-2*exp(-TI/t1) + exp(-TR/t1)).
%
% OPTIONAL INPUTS:
%
% options - algorithm options
%
%
%
%
% Author: Paddy Slator, p.slator@ucl.ac.uk
%
%
% LICENSE
% <inspect toolbox for qMRI analysis>
% Copyright (C) <2020>  <Paddy J. Slator, p.slator@ucl.ac.uk>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

tic;

%% preprocess the image/images

[~,imgfilename,mask,nimg,allimg,imgind,voxind,nvox,nx,ny,nz] = inspect_preprocess_img(img,mask);


%% extract the MR acqusition parameters 

if ischar(gradechoinv)%check if gradechoinv is a path to a file
    gradechoinvfilename = gradechoinv;
    gradechoinv = importdata(gradechoinvfilename);
end


%% unpack algorithm options

%get the default options
default_options = default_options_inspect_map(kernel,gradechoinv,imgfilename); 

if nargin < 5 %if no user defined options    
    options=default_options;
    disp('Using default options.')                
elseif nargin == 5 %amend any user defined options
    options = append_options(options,default_options);
end

%print the options that are going to be used
options
options.ILT

if options.save
    %print save directory, and create it if it doesn't exist
    if exist([options.save_path options.dirname], 'dir')
        disp(['Results will be saved at: ' options.save_path options.dirname])
    else
        mkdir([options.save_path options.dirname])
        disp(['Created directory for saving results at: ' options.save_path options.dirname])
    end
    disp(['Output filenames will end in '  strjoin(options.scan_names,'_')])
end




%% unpack a few well-used options

%number of spectral components to fit
ncomp = options.ncomp;

%turn variables on/off for debugging
onweights = options.onweights;
onF = options.onF;

%the grid on which to estimate the diffusion/relaxometry spectrum
Nk = options.ILT.Nk;
mink = options.ILT.mink;
maxk = options.ILT.maxk;



%% get the kernel dictionary values by doing a dummy fit to the first voxel

%this turns off the actual ILT calculation
options.ILT.onILT = 0; 
ILT_output_test = ILT(allimg(1,:)', gradechoinv, options.ILT);
%store the values for subsequent ILT calls
options.ILT.K = ILT_output_test.K;
options.ILT.grid = ILT_output_test.grid;
if options.ILT.reg
    options.ILT.Kalpha = ILT_output_test.Kalpha;
end
%extract dictionary values required for (non-ILT) calculations
K = ILT_output_test.K;
if options.ILT.reg
    Kalpha = ILT_output_test.Kalpha;
else
    Kalpha = K;
end
%normal fitting for all subsequent calls
options.ILT.onILT=1;




%% initialise the canonical spectra and weights

if strcmp(options.init,'kmeans')
    %intialise clustering with k-means of the signal
    %not sure if this is any good for these - because the signal magnitudes are
    %so different - is there a better clustering metric that we can use?
    %rearrange the image
    roivec = kmeans(allimg,ncomp);

    %convert rois to weights
    weights = zeros([nvox ncomp]);
    for i=1:nvox
        weights(i,roivec(i)) = 1;
    end

    if onF
        %initialise Fcomp values by fitting to weighted signal
        %(not really sure if this is the best way to do it?)
        Fcompvec = cell(1,ncomp);
        Fcomp = cell(1,ncomp);
        for i=1:ncomp
            weightedsig = sum(bsxfun(@times, weights(:,i), allimg));

            %calculate the ILT on the weighted signal for this component
            %ILT_output = ILT_2D(weightedsig, gradechoinv, options.ILT_mean);
            %Fcompvec{i} = ILT_output.Fvec;
            %Fcomp{i} = ILT_output.F;
            
            ILT_output = ILT(weightedsig, gradechoinv, options.ILT_mean);
            Fcompvec{i} = ILT_output.Fvec;
            Fcomp{i} = ILT_output.F;
            
        end

        %calculate spectrum on the overall mean signal
        meansig = mean(allimg);
        ILT_mean = ILT_2D(meansig,gradechoinv,options.ILT_mean);

        %normalise so each spectrum component integrates to the same value
        %as the mean spectrum
        for i=1:ncomp
           Fcomp{i} = sum(ILT_mean.F(:)) * Fcomp{i}./sum(Fcomp{i}(:));
           Fcompvec{i} = sum(ILT_mean.F(:)) * Fcompvec{i}./sum(Fcompvec{i});
        end

        %[Fcompvec, Fcomp, dmean, Cmean] = inspect_map_F_update(allimg,Fcompvec,weights,K,Kalpha,gradechoinv,options);
    else
        Fcomp = options.Fcomp;
        for i=1:length(Fcomp)
            Fcompvec = Fcomp{i}(:);
        end
    end
    disp('initialised with kmeans clustering')
elseif strcmp(options.init,'meanspectrum')
    %initialise spectral components from a fit to the mean signal across
    %the whole image, then calculate weights based on these
    if onF
        %calculate spectrum on the overall mean signal
        meansig = mean(allimg);       
        ILT_mean = ILT(meansig,gradechoinv,options.ILT_mean);       
        output.ILT_mean = ILT_mean;

        %identify the separated peaks in this, and use these as the intial spectra      
        [peaks,peaksvec] = extract_components_from_spectrum(ILT_mean.F);

        %if there are not enough peaks - do something
        if length(peaks) < ncomp

            %randomly combine the peaks to make new peaks
            disp('combining peaks to make new intial compartments')
            while length(peaks) < ncomp
                to_add = randperm(length(peaks),2);

                peaks{end+1} = 0.5 * (peaks{to_add(1)} + peaks{to_add(2)});
                peaksvec{end+1} = 0.5 * (peaksvec{to_add(1)} + peaksvec{to_add(2)});
            end

        end

        %find the ncomp biggest peaks
        for i=1:length(peaks)
            peakweight(i) = sum(peaksvec{i});
        end
        [~, peakweightorder] = sort(peakweight, 'descend');

        %only use these peaks
        Fcomp = peaks(peakweightorder(1:ncomp));
        Fcompvec = peaksvec(peakweightorder(1:ncomp));

        %normalise so each spectrum component integrates to the same value
        %as the mean spectrum
        for i=1:ncomp
           Fcomp{i} = sum(ILT_mean.F(:)) * Fcomp{i}./sum(Fcomp{i}(:));
           Fcompvec{i} = sum(ILT_mean.F(:)) * Fcompvec{i}./sum(Fcompvec{i});
        end

        %reorder the peaks by the first dimension (x-axis)
        maxindcell = cell(ncomp,length(Nk));
        
        for i=1:ncomp
            %find maximum value in this peak
            [~,maxindvec(i)] = max(Fcompvec{i});
                        
            [maxindcell{i,:}] = ind2sub(Nk, maxindvec(i));
                               
        end                   
      
        %order the last dimension
        [~,comporder] = sortrows(maxindcell,length(Nk));
        
        for i=1:ncomp
            Fcomptemp{i} = Fcomp{i};
        end
        for i=1:ncomp
            Fcomp{i} = Fcomptemp{comporder(i)};
            Fcompvec{i} = Fcomp{i}(:);
        end                      
    else
        Fcomp = options.Fcomp;
        for i=1:length(Fcomp)
            Fcompvec{i} = Fcomp{i}(:);
        end
    end
    initvals.Fcomp = Fcomp;
    initvals.Fcompvec = Fcompvec;   
        
    if onweights
        %initialise weights given initial spectra
        currentweights = []; %random starting points               
        weights = inspect_map_weights_update(allimg,currentweights,Fcomp,Kalpha,gradechoinv,options);              
    else
        weights = options.weights;
    end
    disp('initialised spectra and weights by selecting spectral components from fit to mean signal')
elseif strcmp(options.init,'user')

    Fcomp = options.Fcomp;
    for i=1:length(Fcomp)
        Fcompvec{i} = Fcomp{i}(:);
    end

    weights = options.weights;

    initvals.Fcomp = Fcomp;
    initvals.Fcompvec = Fcompvec;

    disp('initialised from user provided values')
elseif strcmp(options.init,'random')
    if onweights
        weights = rand([nvox ncomp]);
        %normalise
        weights = weights./repmat(sum(weights,2),[1 ncomp]);
    else
        weights = options.weights;
    end

    if onF
        %initialise Fcomp values by fitting to weighted signal
        %(not really sure if this is the best way to do it?)
        Fcompvec = cell(1,ncomp);
        Fcomp = cell(1,ncomp);
        for i=1:ncomp
            weightedsig = sum(bsxfun(@times, weights(:,i), allimg));

            %calculate the ILT on the weighted signal for this component
            ILT_output = ILT(weightedsig, gradechoinv, options.ILT_mean);
            
            Fcompvec{i} = ILT_output.Fvec;
            Fcomp{i} = ILT_output.F;
        end

    else
        Fcomp = options.Fcomp;
        for i=1:length(Fcomp)
            Fcompvec = Fcomp{i}(:);
        end
    end

    disp('initialised randomly')
end



%store the initial volume fraction image
for i=1:nimg
    imgweights{i} = voxel_to_image(weights(imgind == i,:),...
        voxind{i},...
        [nx{i} ny{i} nz{i}]);

    initvals.imgweights{i} = imgweights{i};
end

initvals.weights = weights;
initvals.Fcomp = Fcomp;



%preallocate a cell to store the results at each iteration
iter = cell(options.maxiter,1);
if onF && onweights
    for j=1:options.maxiter
        iter{j} = cell(ncomp,1);
    end
end
%preallocate some other things
stepwiselogli = zeros(options.maxiter,1);
stepwiseRESNORM = zeros(options.maxiter,1);

converged = 0;

k = 1;
while k <= options.maxiter && ~converged

    %% main update procedure:
    % 1. update first canonical spectrum,
    % 2. update all voxel weights
    % 3. update second canonical spectrum,
    % 4. update all voxel weights
    % etc. until all spectra have been updated and weights have been updated ncomp times

    if onF
        %spectrum maximisation step for the canonical spectrum components
        for j=1:ncomp %loop through components

            %only update this spectral component
            options.includecomp = zeros(ncomp,1);
            options.includecomp(j) = 1;
            %calculate spectrum for this canonical spectral component,
            %given the other spectral components and weights            
            [Fcompvec, RESNORM, Fcomp, dmean{j}, Cmean{j}] = inspect_map_F_update(allimg,Fcompvec,weights,K,Kalpha,gradechoinv,options);


            if onweights %now update the weights

                %potential speed up - only update voxels where the tolerance is
                %above some limit - at the moment we update all voxels
                %options.updatevox = abs(weightstol(:)) < options.weightstol;

                %do the weight update
                weights = inspect_map_weights_update(allimg,weights,Fcomp,Kalpha,gradechoinv,options);
                              
                %store the parameters at this iteration
                %spectra in two formats
                iter{k}{j}.Fcompvec = Fcompvec;
                iter{k}{j}.Fcomp = Fcomp;
                %weights
                iter{k}{j}.weights = weights;
                %residual    
                iter{k}{j}.Mstep_F_RESNORM = RESNORM;
                %for debugging!
                iter{k}{j}.dmean = dmean;
                iter{k}{j}.Cmean = Cmean;

            else %if debugging/testing just store the user-defined weights
                weights = options.weights;

                iter{k}{j}.Fcompvec = Fcompvec;
                iter{k}{j}.Fcomp = Fcomp;
                iter{k}{j}.Mstep_F_RESNORM = RESNORM;
                iter{k}{j}.weights = weights;

                %for debugging!
                iter{k}{j}.dmean = dmean;
                iter{k}{j}.Cmean = Cmean;
            end

        end
    else
        Fcomp = options.Fcomp;
    end


    if onweights && ~onF %only if we are debugging/testing the weight updates
        %update the weights

        weights = inspect_map_weights_update(allimg,weights,Fcomp,Kalpha,gradechoinv,options);

        %update the final step sizes for the next run
        options.hill.stepsize = stepsize;

        j=1; %only one update per iteration     
        %store the parameters at this iteration
        iter{k}{j}.Fcompvec = Fcompvec;
        iter{k}{j}.Fcomp = Fcomp;
        %store the weights
        iter{k}{j}.weights = weights;
    end

    %go from weights back to image space
    for i=1:nimg
        for j=1:length(iter{k})
            imgweights{i} = voxel_to_image(iter{k}{j}.weights(imgind == i,:),...
                voxind{i},...
                [nx{i} ny{i} nz{i}]);

            iter{k}{j}.imgweights{i} = imgweights{i};
        end
    end


    %calculate the log-likelihood at the end of this iteration
    %augment the signal 
    allimgaug = [allimg zeros(nvox, prod(Nk))];
    
    [logli,RESNORM] = calculate_map_logli_all_voxels(allimg,allimgaug,Fcompvec,Kalpha,weights,gradechoinv,options);

    iter{k}{end}.logli = logli;
    iter{k}{end}.RESNORM = RESNORM;

    stepwiselogli(k) = sum(logli);
    stepwiseRESNORM(k) = RESNORM;


    output.iter = iter;

    disp(['iteration ' num2str(k) ' complete (max iterations ' num2str(options.maxiter) ')'])

    %CONVERGENCE CHECK 
    if k>1
        %check the tolerance of the parameter maps
        weightstol = iter{k-1}{end}.weights - iter{k}{end}.weights;

        iter{k}{end}.weightstol = weightstol;

        if max(abs(weightstol(:))) < options.weightstol
            converged = 1;
            disp(['converged on the ' num2str(k) 'th iteration'])

            conviter = k;

            iter = iter(1:conviter);
            
            stepwiselogli = stepwiselogli(1:conviter);
            stepwiseRESNORM = stepwiseRESNORM(1:conviter);
        end
    end
    %if it went to the maximum iteration
    if k == options.maxiter
        conviter = k;
    end


    k = k+1;
end


runtime=toc;



relabel = options.relabel;

if relabel
    if onF
        %relabel the spectral components based on mean t2* value
        for i=1:options.ncomp
            %get this spectrum
            thisF = iter{conviter}{end}.Fcomp{i};
            %find the mode of the spectrum
            [~,maxindex] = max(thisF(:));
            [~,modepos(i)] = ind2sub(size(thisF), maxindex);
        end
        [~,clusterindex] = sort(modepos);
    end
end

%
if onF       
    %put the spectra into a single mat file
    F = cell(ncomp,1);
    w = cell(ncomp,1);

    if relabel
        for j=1:ncomp
            %put into the relabelled order               
            F{clusterindex(j)} = iter{conviter}{end}.Fcomp{j};
            w{clusterindex(j)} = ILT_output_test.grid;
        end
    else
        for j=1:ncomp
            F{j} = iter{conviter}{end}.Fcomp{j};
            w{j} = ILT_output_test.grid;          
        end
    end
    %make a structure
    spectra.F = F;
    spectra.w = w; 
    output.spectra = spectra;

    if isfield(options,'save') %save the spectra as a mat file
        if options.save
            save([options.save_path options.dirname '/spectra_' num2str(ncomp) '_comp_' strjoin(options.scan_names,'_')],'spectra');
        end
    end
    if relabel %relabel the imgweights
        for i=1:nimg
            imgweightstemp = iter{conviter}{end}.imgweights{i};
            for j=1:ncomp
                imgweightstemp(:,clusterindex(j)) = iter{conviter}{end}.imgweights{i}(:,j);
            end
            imgweights{i} = imgweightstemp;
        end
    end
end

%save the volume fraction image as nifti
if isfield(options,'save')
    if options.save
        for i=1:nimg %separate nifti for each image        
            niftiwrite(imgweights{i},[options.save_path options.dirname '/inspect_map_'  num2str(ncomp) '_comp_' options.scan_names{i} '.nii.gz'])
        end
        output.fullsavepath = [options.save_path options.dirname];
    end
end




%calculate the AIC and BIC
%number of model parameters
%there are ncomp - 1 weight parameters in each voxel - since they sum to 1
k = (ncomp - 1) * nvox + ncomp*( prod(Nk) );
%sample size
n = nvox * size(gradechoinv,1) ;
LogLi = stepwiselogli(end);

AIC = 2*k - 2*LogLi;
BIC = k*log(n) - 2*LogLi;

output.AIC = AIC;
output.BIC = BIC;

output.runtime = runtime;
output.weights = weights;
output.iter = iter;
output.initvals = initvals;

output.stepwiselogli = stepwiselogli;
output.stepwiseRESNORM = stepwiseRESNORM;

output.imgweights = imgweights;

output.options = options;

output.kernel = kernel;

%save a summary of the output
outputsummary.AIC = output.AIC;
outputsummary.BIC = output.BIC;
outputsummary.runtime = output.runtime;
%don't save this - super expensive!
%outputsummary.iter = output.iter;
outputsummary.stepwiselogli = output.stepwiselogli;
outputsummary.stepwiseRESNORM = output.stepwiseRESNORM;
outputsummary.imgweights = output.imgweights;
outputsummary.options = output.options;
outputsummary.initvals = output.initvals;

if isfield(options,'save')
    if options.save
        try
            save([options.save_path options.dirname '/outputsummary_' num2str(ncomp) '_comp_' strjoin(options.scan_names,'_')],'outputsummary');
        catch %if the file is bigger than 2GB - not ideal but still save the output
            save([options.save_path options.dirname '/outputsummary_' num2str(ncomp) '_comp_' strjoin(options.scan_names,'_')],'outputsummary','-v7.3');
        end
    end
end

%do a plot and save it
try    
    plot_inspect_map(output,options)
catch
    warning('plotting function did not work...')
end

end
