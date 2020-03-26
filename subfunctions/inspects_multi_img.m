function output = inspects_multi_img(img,gradechoinv,mask,options)

%inspect quatitative mapping version
%assumes that the signal from each voxel is a weighted sum of contributions 
%from different spectra

if ~iscell(img) %single image as input - just put into a unit cell
    tempimg = cell(1);    
    tempmask = cell(1);
    
    tempimg{1} = img;
    tempmask{1} = mask;
    
    img =  tempimg;
    mask = tempmask;
    
end
    
%multiple images
nimg = length(img);

%have to assume that all images have the same gradechoinv
%not sure if there is an easy way to deal with different ones?

%get the b-values in \times 10^-3 s/mm^2 (nice ADC output)
b = gradechoinv(:,4) * 10^-3;
%get the echo times in seconds
te = gradechoinv(:,5) * 10^-3;

%get the number of diffusion-relaxometry volumes
Nmeas = size(gradechoinv,1);



for i=1:nimg
    %dimensions of each volume
    [nx{i},ny{i},nz{i},~] = size(img{i}(:,:,:,1));
    %number of voxels to fit
    nvox(i)=nnz(mask{i});
end


%concatenate all the images into a single image?
%definately easiest - do it as one big thing "voxel style"

imgvox = cell(nimg,1);
voxind = cell(nimg,1);

allimg = [];
imgind = [];
for i=1:nimg
    %put image into voxel form
    [imgvox{i}, voxind{i}] = image_to_voxel(img{i},mask{i});
    allimg = [allimg; imgvox{i}];
    %store the img index for these voxels
    imgind = [imgind; i*ones(size(imgvox{i},1),1)];
end


%dimensions of allimg
allimgsize = size(allimg);
nvox = allimgsize(1);

tic;


%unpack some other stuff
Nk1 = options.ILT.Nk1;
Nk2 = options.ILT.Nk2;

nclus = options.EM.n_clusters;

%matrix of all possible spectra assignments
allspectcombs = dec2bin(2^nclus-1:-1:0)-'0';
%remove the row with all zeros
allspectcombs = allspectcombs( sum(allspectcombs,2)~=0,:);
%number of possible spectra assigment combinations
%all zeros is not a legit clustering!
ncombs = 2^nclus - 1;
    


if strcmp(options.EM.init,'kmeans')
    %intialise clustering with k-means of the signal
    %not sure if this is any good for these - because the signal magnitudes are
    %so different - is there a good clustering metric that we can use?
    %rearrange the image
    roivec = kmeans(allimg,ncombs);
    
    %convert rois to weights
    weights = zeros([nvox ncombs]);
    for i=1:nvox
        weights(i,roivec(i)) = 1;
    end
    disp('initialised with kmeans clustering')
elseif strcmp(options.EM.init,'meanspectrum')
    %initialise spectral components from a fit to the mean signal across
    %the whole image, then calculate weights based on these
    
    %calculate spectrum on the overall mean signal
    meansig = mean(allimg);
    ILT_mean = ILT_2D(meansig,gradechoinv,options.ILT);   
    output.ILT_mean = ILT_mean;
    Kalpha = ILT_mean.Kalpha;
            
    %identify the separated peaks in this, and use these as the intial spectra        
    binF = imbinarize(ILT_mean.F);        
    connections = bwconncomp(binF);    
    for i=1:connections.NumObjects
        peaks{i} = connections.PixelIdxList{i};
    end        
    %get the starting values of Fvec - the biggest n cluster
    for i=1:nclus
        Fvec{i} = zeros(size(ILT_mean.Fvec));
        Fvec{i}(peaks{i}) = ILT_mean.Fvec(peaks{i});
        
        %[x,y] = ind2sub(peaks{i},size(ILT_mean.F));
        %F{i} = zeros(size(ILT_mean.Fvec));
    end
        
    %initialise with random spectral weights assignments
    p = ones(nvox,nclus);
    for i=1:nvox
        p(i,:) = p(i,:)./sum(p(i,:));
    end
    
    %initialise posterior latent state weights 
    weights = inspects_multi_img_Estep(allimg,gradechoinv,p,Fvec,Kalpha,options);
    
    
    output.initFvec=Fvec;
    output.initp = p;
    output.initweights = weights;
    
    
    disp('initialised spectra and weights by selecting spectral components from fit to mean signal')
else
    %randomly assign initial cluster assignment weights - need to calculate
    %for every possible cluster assigment - there are 2^Nclus-1 of these
    weights = randi(0:1,[nvox ncombs]);
    disp('initialised randomly')
    %add a one in a random row of any all zero rows
    weights(sum(weights,2) == 0, randi(ncombs)) = 1;
    %normalise
    weights = weights./repmat(sum(weights,2),[1 ncombs]);
end





%cluster weights
%for j=1:nclus
%    clusweights(:,j) = sum(weights(:, allspectcombs(:,j)==1 ),2);
%end


%
F = cell(nclus,1);


EMstep = cell(options.EM.n_steps,1);

for k = 1:options.EM.n_steps
    
    %M-step
    %update voxelwise and cluster weights
    p = zeros(nvox,nclus);
    for i=1:nvox         
        for j=1:nclus
            %weight associated with this voxel and spectrum
            %sum over all weights where cluster j = 1
            Nj = sum( weights(i, allspectcombs(:,j)==1 ));                        
            %sum the weights of this cluster across all images       
            p(i,j) = Nj/ncombs ;   
            
             %sum over all weights where cluster j = 1
            %Nj = clusweights(i, j );                        
            %sum the weights of this cluster across all images       
            %p(i,j) = Nj/nclus ;   
            
            
        end
        %normalise within voxel
        p(i,:) = p(i,:)/(sum(p(i,:)));
    end
    
    %calculate spectrum for each cluster based on weights
    ILT_output = cell(nclus,1);
    
    %get the data for each cluster - average signal across all voxels - weighted by
    %the probability that this voxel is in this cluster
    nvol = size(allimg,2);
    
    for j=1:nclus
        %weight the signal in all voxels - only including weights for z_n where
        %this cluster is 1
        
        
        Wn = sum(sum(weights(:,allspectcombs(:,j)==1)));        
        
        weightedsig{j}=[];
        for i=1:ncombs
            if allspectcombs(i,j)              
                weightedsig{j} = [weightedsig{j}; repmat(weights(:,i),[1 nvol]).*(allimg)];                
            end
        end    
        
        weightedmeansig{j} = sum(weightedsig{j},1)/Wn;
                
        %OR
        %weightedsig{j} = repmat(clusweights(:,j),[1 nvol]).*allimg;
        %weightedmeansig{j} = sum(weightedsig{j},1) / sum(clusweights(:,j));
        
    end
    
    
    %fit weighted spectrum and get the error (=log likelihood) for each cluster
    for j=1:nclus                                            
        %calculate the T2-diff spectrum
        ILT_output{j} = ILT_2D(weightedmeansig{j},gradechoinv,options.ILT);
                
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
    
    tic; %slow looped version       
    
    weights = inspects_multi_img_Estep(allimg,gradechoinv,p,Fvec,Kalpha,options);
    
%     for i=1:nvox
%         %get the signal for this voxel
%         S = allimg(i,:)';
%         %get the SNR for this voxel        
%         %SNR = 1/std(S);
%         %SNR = 10000;        
%         sig = std(S(b == 0 & te == min(te)));
%         
%         %augment the signal
%         S = [S; zeros(options.ILT.Nk1 * options.ILT.Nk2, 1)];
%         
%         unnormweight = zeros(ncombs,1);
%         %for j=1:nclus  %loop through clusters
%         for j=1:ncombs %loop through possible clusterings
%             %get this clustering combination
%             clus = allspectcombs(j,:);
%                
%             %calculate the unnormalised expectation for each cluster
%             %combination 
%             %unnormweight(j) = log(p(j)) + sum( log(normpdf(S, Kalpha*Fvec{j}, 1/SNR) ));
%             %normal scale                                   
%             %for m=1:nclus               
%             %    unnormweightelem(m) = clus(m)*p(i,m)*prod(normpdf(S,Kalpha*Fvec{m},sig));
%             %    unnormweightelem(m)
%             %end            
%             %unnormweight(j) = sum(unnormweightelem);
%             
%             %normalise p depending on the cluster arrangment
%             thisp = zeros(size(p(i,:)));
%             clusindex = find(allspectcombs(j,:));
%             thisp( clusindex ) = p(clusindex);            
%             thisp = thisp./sum(thisp);                       
%             
%             %in log-scale
%             for m=1:nclus  
%                 %disp('p(i,m) on this line needs to be normalised depending on what the cluster arrangement is')
%                 %unnormweightelem(m) = log(clus(m)) + log(p(i,m)) -0.5*log(2*pi*sig^2)  - sum((1/(2*sig))*(S - Kalpha*Fvec{m}).^2);                
%                 unnormweightelem(m) = log(clus(m)) + log(thisp(m)) -0.5*log(2*pi*sig^2)  - sum((1/(2*sig))*(S - Kalpha*Fvec{m}).^2);             
%             end
%                                     
%             %allspectcombs(j,:)
%             %unnormweightelem
%             
%             %log-sum-exp trick
%             logsumexp_unnormweightelem = max(unnormweightelem) + log(sum(exp(unnormweightelem - max(unnormweightelem))));            
%                        
%             unnormweight(j) =  logsumexp_unnormweightelem;
%                            
%             %sometimes underflows to NaN if this cluster combination has very low probability
%             if isnan(unnormweight(j))
%                 unnormweight(j) = -Inf;
%             end
%             
%             %unnormweight(j) = log(sum(exp(unnormweightelem)));
%         end
%         %unnormweight
%                              
%         %normalise to get the posterior
%         %do the log-sum-exp trick
%         logsumexp_unnormweight = max(unnormweight) + log(sum(exp(unnormweight - max(unnormweight))));
%         weights(i,:) = exp ( unnormweight - logsumexp_unnormweight );
%         
%        
%                 
%         %get cluster weights
%         for j=1:nclus
%             clusweights(:,j) = sum(weights(:, allspectcombs(:,j)==1 ),2);
%         end
%         
%         %weights(i,:) = unnormweight./sum(unnormweight);
%     end
    toc;
    
   
%     %memory heavy vectorised version!
%     %get the SNR for all voxels        
%     %SNR = (1./std(allimg'))';
%     %sigvec = std(allimg');    
    sigvec = std(allimg(:, b == 0 & te == min(te) )');
    sigvec = sigvec';
%     
%     %augment the signal
%     %augallimg = [allimg zeros(nvox, options.ILT.Nk1 * options.ILT.Nk2)];       
%     %disp('augallimg')
%     %size(augallimg)
%                   
%     %calculate the unnormalised expectation for each cluster
%     unnormweightvec = zeros([nvox nclus]);
%     
%     gaussterms = zeros(nvox,nclus);
%       
%     tic;
%     for j=1:nclus
%         %get the number of elements in the spectra
%         nspectelem = Nk1 * Nk2 + size(allimg,2);
%         
%         
%         %construct matrix of sigma values - each row is the SNR for
%         %that voxel repeated
%         %sigmat =  repmat( sigvec, [1 nspectelem]) ;
%         
%         %disp('sigmat')
%         %size(sigmat)
%         
%         
%         %calculate the gaussian part of the weights in log scale
%         %gaussterms = sum( -1./(2.*sigmat) .* (augallimg - repmat((Kalpha*Fvec{j})',[nvox 1])).^2, 2);
%         
%         %do whole calculation on one messy line to save allocating huge arrays
%         %gaussterms = sum( -1./(2.*repmat( sigvec, [1 nspectelem])) .* ...
%         %    ([allimg zeros(nvox, Nk1 * Nk2)] - ...
%         %    repmat((Kalpha*Fvec{j})',[nvox 1])).^2, 2);
%               
%         if nvox > 10^5 %do the calculation in blocks to save on memory and time
%             %hard code the number of blocks
%             nblocks = 10;
%             blocksize = floor(nvox/nblocks);
%             
%             
%             for i=1:nblocks
%                 if i==nblocks %last block will need to be bigger unless nvox divides 10
%                     block = ((i-1)*blocksize + 1):nvox;
%                     blocksize = length(block);
%                 else
%                     block = ((i-1)*blocksize + 1) : i*blocksize;
%                 end
%                 
%                 %do whole calculation on one messy line to save allocating huge arrays
%                 gaussterms(block,j) =...
%                     sum( -1./(2.*repmat( sigvec(block).^2, [1 nspectelem])) .* ...
%                     ([allimg(block,:) zeros(blocksize, Nk1 * Nk2)] - ...
%                     repmat((Kalpha*Fvec{j})',[blocksize 1])).^2, 2);
%             end
%         else
%             %do whole calculation on one messy line to save allocating huge arrays
%                 gaussterms(:,j) =...
%                     sum( -1./(2.*repmat( sigvec.^2, [1 nspectelem])) .* ...
%                     ([allimg zeros(nvox, Nk1 * Nk2)] - ...
%                     repmat((Kalpha*Fvec{j})',[nvox 1])).^2, 2);
%         end
%         
%         
%         %calculate the unnormalised expectation for each cluster in log scale
%         unnormweightvec(:,j) = log(p(j)) + gaussterms(:,j);
%         
%     end
%     toc; 
%     
%     
%     %normalise to get the posterior 
%     %do the log-sum-exp trick
%     %get max of the weights for each voxel
%     maxweights = max(unnormweightvec,[],2);    
%     logsumexp_unnormweightvec = max(unnormweightvec,[],2) + log(sum(exp(unnormweightvec -  repmat(maxweights,[1 nclus]) ), 2));
%                   
%     weightsvec = exp( unnormweightvec - repmat(logsumexp_unnormweightvec, [1 nclus]) );
%     
%     %just in case want to test against voxelwise code
%     weights = weightsvec;
                   
    
    
    EMstep{k}.weights = weights;
    %EMstep{k}.clusweights = clusweights;
    EMstep{k}.p = p;
    EMstep{k}.ILT_output = ILT_output;
 
    
    %go from weight probabilities to a maximum likelihood roi for each
    %image
    for i=1:nimg                
        imgweights{i} = voxel_to_image(weights(imgind == i,:),...
            voxind{i},...
            [nx{i} ny{i} nz{i}]);
        
        MLroi{i} = weights_to_roi(imgweights{i},mask{i});
        
        EMstep{k}.imgweights{i} = imgweights{i};
        EMstep{k}.MLroi{i} = MLroi{i};        
    end
    
    %calculate the log-likelihood
    logliterms=zeros(nvox,nclus);
    for i=1:nvox
        S = allimg(i,:)';
        %get the SNR for this voxel                     
        sig = std(S(b == 0 & te == min(te)));
        %augment the signal
        S = [S; zeros(options.ILT.Nk1 * options.ILT.Nk2, 1)];
        for j=1:nclus                                               
            logliterms(i,j)  = log(p(i,j)) -0.5*log(2*pi*sig^2)  - sum((1/(2*sig))*(S - Kalpha*Fvec{j}).^2);
        end
    end
    
    %for j=1:nclus
    %    logliterms(:,j) = log(p(j)) -0.5*log(2*pi*sigvec.^2) + gaussterms(:,j);
    %    %logliterms(:,j) = log(p(j)) - gaussterms(:,j)
    %end
    
    %do the log-sum-exp trick
    %get max of the weights for each voxel
    maxlogliterms = max(logliterms,[],2);    
    logli = sum ( max(logliterms,[],2) + log(sum(exp(logliterms -  repmat(maxlogliterms,[1 nclus]) ), 2)));
    
    
    EMstep{k}.logli = logli;
    stepwiselogli(k) = logli;
    
    disp(['EM step ' num2str(k) ' of ' num2str(options.EM.n_steps) ' complete'])     
    
    
    
end


runtime=toc;







%relabel the clusters based on mean t2* value
for i=1:options.EM.n_clusters
    %get this spectrum
    thisF = EMstep{end}.ILT_output{i}.F;
    %find the mode of the spectrum
    [~,maxindex] = max(thisF(:));
    [~,modepos(i)] = ind2sub(size(thisF), maxindex);
end
[~,clusterindex] = sort(modepos);

if isfield(options,'save')
    if options.save
        %save everything in its own directory
        dirname = [num2str(options.EM.n_clusters) '_clusters_' ...
            num2str(options.EM.n_steps) '_steps_' ...
            'alpha_' num2str(options.ILT.alpha)];
        
        mkdir([options.save_path dirname])
    end
end

%save the spectra in a single mat file
F = cell(options.EM.n_clusters,1);
w1 = cell(options.EM.n_clusters,1);
w2 = cell(options.EM.n_clusters,1);
for j=1:options.EM.n_clusters
    %save in the relabelled order
    F{clusterindex(j)} = EMstep{end}.ILT_output{j}.F;
    w1{clusterindex(j)} = EMstep{end}.ILT_output{j}.w1;
    w2{clusterindex(j)} = EMstep{end}.ILT_output{j}.w2;
end
if isfield(options,'save')
    if options.save
        save([options.save_path dirname '/spectra'],'F','w1','w2');
    end
end

%relabel the MLrois
for i=1:nimg
    MLroitemp = MLroi{i};
    for j=1:options.EM.n_clusters
        MLroitemp(MLroi{i} == clusterindex(j)) = j;
    end
    MLroi{i} = MLroitemp;
end

%relabel the imgweights
for i=1:nimg
    imgweightstemp = imgweights{i};
    for j=1:options.EM.n_clusters
        imgweightstemp(:,clusterindex(j)) = imgweights{i}(:,j);
    end
    imgweights{i} = imgweightstemp;
end

if isfield(options,'save')
    if options.save
        %save the ROIs as nifti files
        %
        for i=1:nimg
            roinii = make_nii(MLroi{i});
            save_nii(roinii,[options.save_path dirname '/MLroi_' options.scan_names{i} '.nii.gz'])
        end
        
        output.fullsavepath = [options.save_path dirname];
    end
end

%now relabel the clusters in everything for the output
tempweights = weights;
for j=1:options.EM.n_clusters
    tempweights(:,j) = weights(:,clusterindex(j));
end
weights = tempweights;

tempp = p;
for j=1:options.EM.n_clusters
    tempp(j) = p(clusterindex(j));
end
p = tempp;


for k = 1:options.EM.n_steps
    tempILT_output = EMstep{k}.ILT_output;
    for j=1:options.EM.n_clusters
        tempILT_output{j} = EMstep{k}.ILT_output{clusterindex(j)};
    end
    EMstep{k}.ILT_output = tempILT_output;
end
%the final EM step
ILT_output = tempILT_output;
        


%calculate the AIC and BIC
%number of model parameters
%k = 2*nclus;
k = 2 * (nclus + nclus*(Nk1*Nk2));
%sample size
n = nvox * Nmeas;
LogLi = stepwiselogli(end);

AIC = 2*k - 2*LogLi;
BIC = k*log(n) - 2*LogLi;

output.AIC = AIC;
output.BIC = BIC;
        

output.runtime = runtime;
output.weights = weights;
output.p = p;
output.ILT_output = ILT_output;
output.EMstep = EMstep;
%save the model fit object in the pwd


output.stepwiselogli = stepwiselogli;

output.MLroi = MLroi;
output.imgweights = imgweights;




output.allspectcombs = allspectcombs;

end









