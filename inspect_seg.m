function [output,outputsummary] = inspect_seg(img,gradechoinv,mask,kernel,options)

%inspect segmentation version
%assumes that each voxel is associated with a single cluster, and each
%cluster has an associated spectrum
%
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

%% preprocess the image/images

[~,imgfilename,mask,nimg,allimg,imgind,voxind,nvox,nx,ny,nz] = inspect_preprocess_img(img,mask);


%% unpack algorithm options
 
%get the default options
default_options = default_options_inspect_seg(kernel,gradechoinv,imgfilename); 

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


%% extract the MR acqusition parameters 

if ischar(gradechoinv)%check if gradechoinv is a path to a file
    gradechoinvfilename = gradechoinv;
    gradechoinv = importdata(gradechoinvfilename);
end




%% get the kernel dictionary values by doing a dummy fit to the first voxel
disp('Calculating the kernel dictionary for all ILT fits.')
tic;
%this turns off the actual ILT calculation
options.ILT.onILT = 0; 
%do the dummy ILT calculation
ILT_output_test = ILT(allimg(1,:)', gradechoinv, options.ILT);
%store the values
options.ILT.K = ILT_output_test.K;
options.ILT.grid = ILT_output_test.grid;
if options.ILT.reg
    options.ILT.Kalpha = ILT_output_test.Kalpha;
end
%normal fitting for all subsequent calls
options.ILT.onILT=1;
time=toc;
disp(['Calculated the kernel dictionary for all ILT fits. It took ' num2str(time) ' seconds.'])




%% unpack a few well-used options

%grid size
Nk = options.ILT.Nk;

%min/max grid values
mink = options.ILT.mink;
maxk = options.ILT.maxk;

%number of clusters
nclus = options.nclus;

%number of image volumes
Nmeas = size(gradechoinv,1);




%% initialise spectra and weights

if strcmp(options.init,'kmeans')
        %intialise clustering with k-means of the signal
        %not sure if this is any good for these - because the signal magnitudes are
        %so different - is there a good clustering metric that we can use?
        %rearrange the image
        disp('Initialising the voxelwise clusters with kmeans.') 
        tic;
        roivec = kmeans(allimg,nclus);        
        %convert rois to weights
        weights = zeros([nvox nclus]);
        for i=1:nvox
            weights(i,roivec(i)) = 1;
        end                
        initvals.weights = weights;
        time=toc;
        disp(['Initialised the voxelwise clusters with kmeans. It took ' num2str(time) ' seconds.']) 
else
    %randomly assign initial cluster weights - final dimension is
    %weights for each cluster
    
    %randomly 0 or 1
    weights = randi(0:1,[nvox nclus]);    
    %prevent divide 0 errors    
    weights = weights + eps;
    %normalise
    weights = weights./repmat(sum(weights,2),[1 nclus]);
    disp('initialised randomly (0 or 1)')
    
%     %random in [0,1]
%     weights = rand([nvox nclus]);
%     %normalise
%     weights = weights./repmat(sum(weights,2),[1 nclus]);
%     disp('initialised randomly in [0,1]')
   
        
        
    %weights = (1/nclus)*ones([nvox nclus]);
    %disp('initialised weights uniformly across all clusters!')

    initvals.weights = weights;
end



%
F = cell(nclus,1);

%% main EM routine

iter = cell(options.maxiter,1);
converged = 0;

k = 1;
while k <= options.maxiter && ~converged
        
    %M-step
    %update cluster probabilities
    
    for j=1:nclus
        %sum the weights of this cluster across all images       
        p(j) = (1/nvox) *  sum(weights(:,j));                       
    end
    %normalise
    p = p/(sum(p));
    
    
    %calculate spectrum for each cluster based on weights
    ILT_output = cell(nclus,1);
    
    %get the data for each cluster - average signal across all voxels - weighted by
    %the probability that this voxel is in this cluster
    nvol = size(allimg,2);
    for j=1:nclus
        %weight the signal in all voxels
        %weightedsig{j} = repmat(weights(:,j),[1 nvol]).*allimg;                
        %normalise
        %weightedmeansig{j} = sum(weightedsig{j},1) / sum(weights(:,j));       
        
        %do this all in one line to save memory          
        weightedmeansig{j} = sum(repmat(weights(:,j),[1 nvol]).*allimg,1) / sum(weights(:,j));
    end
    
    %fit weighted spectrum and get the error (=log likelihood) for each cluster
    for j=1:nclus                                            
        %calculate the T2-diff spectrum
        %ILT_output{j} = ILT_2D(weightedmeansig{j},gradechoinv,options.ILT);
        %arbitrary kernal and dimension       
        disp(['Calculating ILT: EM step ' num2str(k) ' of ' num2str(options.maxiter)...
            ' max, cluster ' num2str(j) '.' ])
        tic;
        ILT_output{j} = ILT(weightedmeansig{j},gradechoinv,options.ILT);
        time=toc;
        disp(['Calculated ILT: EM step ' num2str(k) ' of ' num2str(options.maxiter)...
            ' max, cluster ' num2str(j) '. It took ' num2str(time) ' seconds.'])
        
        %get the error - sum of the squared residuals 
        SSR(j) = ILT_output{j}.RESNORM;
        
        %get the spectrum - vector format
        Fvec{j} = ILT_output{j}.Fvec;
        %matrix
        F{j} = ILT_output{j}.F;
        %this is always the same but just store it here for convenience       
        Kalpha = ILT_output{j}.Kalpha;        
    end
    
    
    %E-step
        
    loop=0;
    
    if loop %slow looped version
        tic; 
        for i=1:nvox
            %get the signal for this voxel
            S = allimg(i,:)';
            %get the SNR for this voxel
            %SNR = 1/std(S);
            %SNR = 10000;
            sig = std(S(b == 0 & te == min(te)));
            
            %augment the signal
            S = [S; zeros(prod(Nk), 1)];            
            
            unnormweight = zeros(nclus,1);
            for j=1:nclus  %loop through clusters
                %calculate the unnormalised expectation for each cluster
                %in log-scale
                %unnormweight(j) = log(p(j)) + sum( log(normpdf(S, Kalpha*Fvec{j}, 1/SNR) ));
                unnormweight(j) = log(p(j)) + sum( - (1/(2*sig))*(S - Kalpha*Fvec{j}).^2);
            end
                        
            
            %normalise to get the posterior
            %do the log-sum-exp trick
            logsumexp_unnormweight = max(unnormweight) + log(sum(exp(unnormweight - max(unnormweight))));
            
            weights(i,:) = exp ( unnormweight - logsumexp_unnormweight );
        end
        toc;
    else   
        %memory heavy vectorised version!
        %get the SNR for all voxels
        %SNR = (1./std(allimg'))';
        %sigvec = std(allimg');      
                     
        %sigvec = std(allimg(:, b == 0 & te == min(te) )');
        %sigvec = sigvec';      
        
        sigvec = estimate_sd(allimg,gradechoinv);       
         
        %augment the signal
        %augallimg = [allimg zeros(nvox, options.ILT.Nk1 * options.ILT.Nk2)];
        %disp('augallimg')
        %size(augallimg)
        
        %calculate the unnormalised expectation for each cluster
        unnormweightvec = zeros([nvox nclus]);
        
        gaussterms = zeros(nvox,nclus);
        
        for j=1:nclus
            %get the number of elements in the spectra        
            %nspectelem = prod(Nk) + size(allimg,2);
            
            
            %construct matrix of sigma values - each row is the SNR for
            %that voxel repeated
            %sigmat =  repmat( sigvec, [1 nspectelem]) ;
            
            %disp('sigmat')
            %size(sigmat)
            
            
            %calculate the gaussian part of the weights in log scale
            %gaussterms = sum( -1./(2.*sigmat) .* (augallimg - repmat((Kalpha*Fvec{j})',[nvox 1])).^2, 2);
            
            %do whole calculation on one messy line to save allocating huge arrays
            %gaussterms = sum( -1./(2.*repmat( sigvec, [1 nspectelem])) .* ...
            %    ([allimg zeros(nvox, Nk1 * Nk2)] - ...
            %    repmat((Kalpha*Fvec{j})',[nvox 1])).^2, 2);
            
            if nvox > 10^5 %do the calculation in blocks to save on memory and time
                %hard code the number of blocks
                nblocks = 10;
                blocksize = floor(nvox/nblocks);
                
                
                for i=1:nblocks
                    if i==nblocks %last block will need to be bigger unless nvox divides 10
                        block = ((i-1)*blocksize + 1):nvox;
                        blocksize = length(block);
                    else
                        block = ((i-1)*blocksize + 1) : i*blocksize;
                    end
                    
                    %do whole calculation on one messy line to save allocating huge arrays
                    if options.ILT.reg
                        nterms = prod(Nk) + size(allimg,2);

                        gaussterms(block,j) =...
                            sum( -1./(2.*repmat( sigvec(block).^2, [1 nterms])) .* ...
                            ([allimg(block,:) zeros(blocksize, prod(Nk))] - ...
                            repmat((Kalpha*Fvec{j})',[blocksize 1])).^2, 2);
                    else
                        nterms = size(allimg,2);
                        
                        gaussterms(block,j) =...
                            sum( -1./(2.*repmat( sigvec(block).^2, [1 nterms])) .* ...
                            (allimg(block,:) - ...
                            repmat((Kalpha*Fvec{j})',[blocksize 1])).^2, 2);                        
                    end
                    
                    
                    
                end
            else
                %do whole calculation on one messy line to save allocating huge arrays               
                if options.ILT.reg
                    %need to augment the data, i.e. [allimg zeros(nvox, prod(Nk))]
                    % 
                    nterms = prod(Nk) + size(allimg,2);
                     
                    gaussterms(:,j) =...
                        sum( -1./(2.*repmat( sigvec.^2, [1 nterms])) .* ...
                        ([allimg zeros(nvox, prod(Nk))] - ...
                        repmat((Kalpha*Fvec{j})',[nvox 1])).^2, 2);
                else                    
                    nterms = size(allimg,2);
                    
                    gaussterms(:,j) =...
                        sum( -1./(2.*repmat( sigvec.^2, [1 nterms])) .* ...
                        (allimg  - ...
                        repmat((Kalpha*Fvec{j})',[nvox 1])).^2, 2);
                end
            end
            
            
            %calculate the unnormalised expectation for each cluster in log scale
            unnormweightvec(:,j) = log(p(j)) + gaussterms(:,j);
            
        end
        
        
        
        %normalise to get the posterior
        %do the log-sum-exp trick
        %get max of the weights for each voxel
        maxweights = max(unnormweightvec,[],2);
        logsumexp_unnormweightvec = max(unnormweightvec,[],2) + log(sum(exp(unnormweightvec -  repmat(maxweights,[1 nclus]) ), 2));
        
        weightsvec = exp( unnormweightvec - repmat(logsumexp_unnormweightvec, [1 nclus]) );
        
        %just in case want to test against voxelwise code
        weights = weightsvec;
    end             
    
    
    iter{k}.weights = weights;
    iter{k}.p = p;
    iter{k}.ILT_output = ILT_output;
 
    
    %go from weight probabilities to a maximum likelihood roi for each
    %image
    for i=1:nimg                
        imgweights{i} = voxel_to_image(weights(imgind == i,:),...
            voxind{i},...
            [nx{i} ny{i} nz{i}]);
        
        MLroi{i} = weights_to_roi(imgweights{i},mask{i});
        
        iter{k}.imgweights{i} = imgweights{i};
        iter{k}.MLroi{i} = MLroi{i};        
    end
    
    if ~loop
        %calculate the log-likelihood
        for j=1:nclus
            logliterms(:,j) = log(p(j)) -0.5*log(2*pi*sigvec.^2) + gaussterms(:,j);
            %logliterms(:,j) = log(p(j)) - gaussterms(:,j)
        end
        
        %normalise to get the posterior
        %do the log-sum-exp trick
        %get max of the weights for each voxel
        maxlogliterms = max(logliterms,[],2);
        logli = sum ( max(logliterms,[],2) + log(sum(exp(logliterms -  repmat(maxlogliterms,[1 nclus]) ), 2)));
                
        iter{k}.logli = logli;
        stepwiselogli(k) = logli;        
    end
    
    
    disp(['EM step ' num2str(k) ' of ' num2str(options.maxiter)  ' max complete'])     
    
    
    %CONVERGENCE CHECK 
    if k>1
        %check the tolerance of the maps
        weightstol = iter{k-1}.weights - iter{k}.weights;

        iter{k}.weightstol = weightstol;

        if max(abs(weightstol(:))) < options.weightstol
            converged = 1;
            disp(['converged on the ' num2str(k) 'th iteration'])

            conviter = k;

            iter = iter(1:conviter);
        end
    end
    %if it went to the maximum iteration
    if k == options.maxiter
        conviter = k;
    end


    k = k+1;
    
    
    
end


runtime=toc;





%% tidy up output

%relabel the clusters based on mean t2* value
for i=1:options.nclus
    %get this spectrum
    thisF = iter{end}.ILT_output{i}.F;
    %find the mode of the spectrum
    [~,maxindex] = max(thisF(:));
    [~,modepos(i)] = ind2sub(size(thisF), maxindex);
end
[~,clusterindex] = sort(modepos);



%save the spectra in a single mat file
F = cell(options.nclus,1);
w = cell(options.nclus,1);

%w1 = cell(options.nclus,1);
%w2 = cell(options.nclus,1);
for j=1:options.nclus
    %save in the relabelled order   
    F{j} = iter{end}.ILT_output{clusterindex(j)}.F;
    w{j} = iter{end}.ILT_output{clusterindex(j)}.grid;
end
if isfield(options,'save')
    if options.save
        save([options.save_path options.dirname '/spectra_' strjoin(options.scan_names,'_')],'F','w');
    end
end

%relabel the MLrois
for i=1:nimg
    MLroitemp = MLroi{i};
    for j=1:options.nclus
        MLroitemp(MLroi{i} == clusterindex(j)) = j;
    end
    MLroi{i} = MLroitemp;
end

%relabel the imgweights
for i=1:nimg
    imgweightstemp = imgweights{i};
    for j=1:options.nclus
        imgweightstemp(:,clusterindex(j)) = imgweights{i}(:,j);
    end
    imgweights{i} = imgweightstemp;
end



%now relabel the clusters in everything for the output
tempweights = weights;
for j=1:options.nclus
    tempweights(:,j) = weights(:,clusterindex(j));
end
weights = tempweights;

tempp = p;
for j=1:options.nclus
    tempp(j) = p(clusterindex(j));
end
p = tempp;


for k = 1:conviter
    tempILT_output = iter{k}.ILT_output;
    for j=1:options.nclus
        tempILT_output{j} = iter{k}.ILT_output{clusterindex(j)};
    end
    iter{k}.ILT_output = tempILT_output;
end
%the final EM step
ILT_output = tempILT_output;
        

if ~loop
    %calculate the AIC and BIC
    %number of model parameters
    %k = 2*nclus;
    k = 2 * (nclus + nclus*(prod(Nk)));
    %sample size
    n = nvox * Nmeas;
    LogLi = stepwiselogli(end);
    
    AIC = 2*k - 2*LogLi;
    BIC = k*log(n) - 2*LogLi;
    
    output.AIC = AIC;
    output.BIC = BIC;
    
    output.stepwiselogli = stepwiselogli;
end

output.runtime = runtime;
output.weights = weights;
output.p = p;
output.ILT_output = ILT_output;
output.iter = iter;
%save the model fit object in the pwd



output.MLroi = MLroi;
output.imgweights = imgweights;

output.initvals = initvals;

output.options = options;

output.kernel = options.kernel;
output.params = GetKernelParameterStrings(options.kernel);



if isfield(options,'save')
    if options.save
        %save the ROIs and img weights as nifti files
        %
        if ~isfield(options,'scan_names') %if no scan name provided leave empty
            options.scan_names = cell(nimg,1);
        end
        for i=1:nimg            
            niftiwrite(MLroi{i},[options.save_path options.dirname '/MLroi_' options.scan_names{i}])
            niftiwrite(imgweights{i},[options.save_path options.dirname '/imgweights_' options.scan_names{i}])
        end
        
        output.fullsavepath = [options.save_path options.dirname];
    end
end


%save a summary of the output - loglikelihood, BIC, etc.
if ~loop
    outputsummary.stepwiselogli = output.stepwiselogli;
    outputsummary.BIC = BIC;
    outputsummary.AIC = AIC;
end

outputsummary.F = F;
outputsummary.w = w;


%outputsummary.weights = weights;
outputsummary.p = p;

outputsummary.initvals = initvals;

outputsummary.options = options;


if isfield(options,'save')
    if options.save
        save([options.save_path options.dirname '/outputsummary_' options.scan_names{i}],'outputsummary');
    end        
end



end









