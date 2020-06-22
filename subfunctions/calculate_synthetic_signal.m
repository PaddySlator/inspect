function [KalphaFvec] = calculate_synthetic_signal(Fcompvec,weights,Kalpha,options)


if iscell(Fcompvec)
    Fcompvec = sparse(cell2mat(Fcompvec));
end

nvox = size(weights,1);
Nk = options.ILT.Nk;

%Fvec = zeros(nvox,Nk);
Fvec = zeros(nvox, prod(Nk));



for i=1:nvox           
    
    %Fvec(i,:) = sum(repmat(weights(i,:),[Nk 1]).*Fcompvec, 2);    
    Fvec(i,:) = sum(repmat(weights(i,:),[prod(Nk) 1]).*Fcompvec, 2);
    
    %Fvecloop(i,:) = sum(bsxfun(@times,weights(i,:), Fcompvec), 2);    
end

%construct matrix where every row is K*Fvec(z_n) for voxel n
%i.e. the synthetic signal
KalphaFvec = (sparse(Kalpha) * sparse(Fvec)')';  





