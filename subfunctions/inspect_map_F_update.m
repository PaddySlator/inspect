function [Fcompvec, RESNORM, Fcomp, dmean, Cmean] = inspect_map_F_update(allimg,FcompvecPrevStep,weights,K,Kalpha,options)

ncomp = options.ncomp;
nvox = size(allimg,1);
Nmeas = size(allimg,2);


%Nk1 = options.ILT.Nk1;
%Nk2 = options.ILT.Nk2;

%Nk = Nk1*Nk2;

Nk = options.ILT.Nk;


%boolean denoting whether to update this compartment
includecomp = options.includecomp;



%calculate the signal components accounted by for by each spectrum


%do this part ignoring regularisation
noreg = 1;
if noreg
    KFvecProd = zeros(Nmeas, nvox,ncomp);
    for j=1:ncomp
        KFmat = repmat(Kalpha * FcompvecPrevStep{j}, [1 nvox] );
        
        weightsmat = repmat(weights(:,j), [1 Nmeas])';         
        
        KFvecProd(:,:,j) = weightsmat .* KFmat;
    end
else
    %do this part with regularisation
    %KalphaFvecProd = zeros(Nmeas + Nk, nvox,ncomp);
    KalphaFvecProd = cell(ncomp,1);
    
    for j=1:ncomp               
        KalphaFmat = repmat(Kalpha * FcompvecPrevStep{j}, [1 nvox] );
        
        %weightsmat = repmat(weights(:,j), [1 Nmeas + Nk])';
        weightsmat = repmat(weights(:,j), [1 Nmeas + prod(Nk)])';
        
        KalphaFvecProd{j} = weightsmat .* KalphaFmat;
    end
end

Fcompvec = cell(1,ncomp);
RESNORM = zeros(1,ncomp);
Fcomp = cell(1,ncomp);

for j=1:ncomp
    if includecomp(j)
        %fast vectorised code for calculating Cterms and dterms
        %indicies of the other components
        othercompindex = [1:(j-1) (j+1):ncomp];
        
        %Cmat = zeros(Nmeas,Nk1 * Nk2, nvox);
        %for i=1:nvox
        %    Cmat(:,:,i) = weights(i,j) * K;
        %end
        %Cmeanfast = sum(Cmat,3);
        
        if noreg %without reg
            Cmeanfast = sum(weights(:,j)) * K;
        else %with reg
            Cmeanfast = sum(weights(:,j)) * Kalpha;
        end
        
        %dterms
        if noreg %without reg
            KFvecProdAll = sum(KFvecProd(:,:,othercompindex),3);
        else %with reg
            KalphaFvecProdAll = sum(KalphaFvecProd{j},3);
        end
        
        if noreg %without reg
            SminusKF = allimg' - KFvecProdAll;
        else %with reg (need to augment the signal)
            %allimgaug = [allimg zeros(nvox, Nk)];
            allimgaug = [allimg zeros(nvox, prod(Nk))];
            SminusKalphaF = allimgaug' - KalphaFvecProdAll;
        end
        
        
        %sum along voxels
        if noreg %without reg
            dmeanfast = sum(SminusKF,2);
        else %with reg
            dmeanfast = sum(SminusKalphaF,2);
        end
        
        %     figure;
        %     plot(allimg')
        %     figure;
        %     plot(KFvecProdAll)
        %     figure;
        %     plot(dmeanfast,'o')
        %     return
        
        dmean = dmeanfast;
        Cmean = Cmeanfast;
        
        if noreg
            if options.ILT.reg
                %augment the dictionary/kernel
                
                %H = speye(Nk1*Nk2,Nk1*Nk2);
                H = speye(prod(Nk), prod(Nk));
                
                %C = [C; options.ILT.alpha*H];
                Cmean = [Cmean; options.ILT.alpha*H];
                %augment the signal
                %d = [d; zeros(Nk1 * Nk2, 1)];
                
                %dmean = [dmean; zeros(Nk1 * Nk2, 1)];                
                dmean = [dmean; zeros(prod(Nk), 1)]; 
            end
        end
        
       
        [Fcompvec{j},RESNORM(j)] = lsqnonneg(Cmean,dmean);
      
        %Fcomp{j} = reshape(Fcompvec{j},[Nk1 Nk2]);       
        
    else
        %don't update this compartment
        Fcompvec{j} = FcompvecPrevStep{j};
        
        %Fcomp{j} = reshape(Fcompvec{j},[Nk1 Nk2]);                
    end
    
    if length(Nk) > 1 %rearrange the spectra if bigger than 1D
        Fcomp{j} = reshape(Fcompvec{j}, Nk);
    else
        Fcomp{j} = Fcompvec{j};
    end
    
end


        
        


end