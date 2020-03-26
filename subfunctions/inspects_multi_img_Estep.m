function weights = inspects_multi_img_Estep(allimg,gradechoinv,p,Fvec,Kalpha,options)

nclus = options.EM.n_clusters;
nvox = size(allimg,1);

%matrix of all possible spectra assignments
allspectcombs = dec2bin(2^nclus-1:-1:0)-'0';
%remove the row with all zeros
allspectcombs = allspectcombs( sum(allspectcombs,2)~=0,:);
%number of possible spectra assigment combinations
%all zeros is not a legit clustering!
ncombs = 2^nclus - 1;


b = gradechoinv(:,4);
te = gradechoinv(:,4);


weights = zeros(nvox,ncombs);


for i=1:nvox;    
    %get the signal for this voxel
    S = allimg(i,:)';
    %get the SNR for this voxel
    %SNR = 1/std(S);
    %SNR = 10000;
    sig = std(S(b == 0 & te == min(te)));
    
    %augment the signal
    S = [S; zeros(options.ILT.Nk1 * options.ILT.Nk2, 1)];
    
    unnormweight = zeros(ncombs,1);
    %for j=1:nclus  %loop through clusters
    for j=1:ncombs %loop through possible clusterings
        %get this clustering combination
        clus = allspectcombs(j,:);
        
        %calculate the unnormalised expectation for each cluster
        %combination
        %unnormweight(j) = log(p(j)) + sum( log(normpdf(S, Kalpha*Fvec{j}, 1/SNR) ));
        %normal scale
        %for m=1:nclus
        %    unnormweightelem(m) = clus(m)*p(i,m)*prod(normpdf(S,Kalpha*Fvec{m},sig));
        %    unnormweightelem(m)
        %end
        %unnormweight(j) = sum(unnormweightelem);
        
        %normalise p depending on the cluster arrangment
        thisp = zeros(size(p(i,:)));
        clusindex = find(allspectcombs(j,:));
        thisp( clusindex ) = p(clusindex);
        thisp = thisp./sum(thisp);
        
        %in log-scale
        for m=1:nclus
            %disp('p(i,m) on this line needs to be normalised depending on what the cluster arrangement is')
            %unnormweightelem(m) = log(clus(m)) + log(p(i,m)) -0.5*log(2*pi*sig^2)  - sum((1/(2*sig))*(S - Kalpha*Fvec{m}).^2);
            unnormweightelem(m) = log(clus(m)) + log(thisp(m)) -0.5*log(2*pi*sig^2)  - sum((1/(2*sig))*(S - Kalpha*Fvec{m}).^2);
        end                    
        
        %clus
        %unnormweightelem
        
        
        %log-sum-exp trick
        logsumexp_unnormweightelem = max(unnormweightelem) + log(sum(exp(unnormweightelem - max(unnormweightelem))));
        
        unnormweight(j) =  logsumexp_unnormweightelem;
                
        
        %sometimes underflows to NaN if this cluster combination has very low probability
        if isnan(unnormweight(j))
            unnormweight(j) = -Inf;
        end
        
        %unnormweight(j) = log(sum(exp(unnormweightelem)));
    end
    %unnormweight'   
    
    %normalise to get the posterior
    %do the log-sum-exp trick
    %logsumexp_unnormweight = max(unnormweight) + log(sum(exp(unnormweight - max(unnormweight))));
    % weights(i,:) = exp ( unnormweight - logsumexp_unnormweight );

    weights(i,:) = exp(unnormweight')./sum(exp(unnormweight));
    
       
    %weights(i,:)
    %get cluster weights
    %for j=1:nclus
    %    clusweights(:,j) = sum(weights(:, allspectcombs(:,j)==1 ),2);
    %end
    
    
    
end