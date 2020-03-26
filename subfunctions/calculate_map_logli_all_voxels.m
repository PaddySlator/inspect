function [LogLi, RESNORM] = calculate_map_logli_all_voxels(allimg,allimgaug,Fcompvec,Kalpha,weights,gradechoinv,options)

%calculate the complete data log-likelihood

%Fcompvec is NK by ncomp 
%weights is nvox by ncomp

%Nk1 = options.ILT.Nk1;
%Nk2 = options.ILT.Nk2;
%Nk = Nk1 * Nk2;
Nk = options.ILT.Nk;


%get the SNR for all voxels
%sigvec = std(allimg(:, b == 0 & te == min(te) )');
%sigvec = sigvec';

sigvec = estimate_sd(allimg,gradechoinv);       


nvox = size(allimg,1);
ncomp = options.ncomp;

% tic;
% %augment the signal
% allimgaug = [allimg zeros(nvox, Nk1 * Nk2)];
% toc;
    



%calculate F(z_n) - the effective voxelwise spectrum - for all voxels at once
%MUCH SLOWER FOR BIG IMAGES     
% tic;
% 
% %convert cell to array
% if iscell(Fcompvec)
%     Fcompvec = cell2mat(Fcompvec);
% end
% 
% %put the Fcomp vectors into nvox by Nk by Ncomp array
% Fcompvecarray = repmat(Fcompvec, [1 1 nvox]);
% Fcompvecarray = permute(Fcompvecarray, [3 1 2]);    
% 
% %put weights into nvox by Nk by Ncomp array
% weights_array = repmat(weights, [1 1 Nk]);
% weights_array = permute(weights_array,[1 3 2]);
% 
% %multiply weights by Fcomp then sum along compartments
% %to construct the effective spectrum for each voxel
% %(nvox by Nk array)
% Fvec = sum(weights_array.*Fcompvecarray,3);
% toc;

%construct F(z_n) - the effective voxelwise spectrum - by looping sparse matrices
%convert cell to array

if iscell(Fcompvec)
    Fcompvec = sparse(cell2mat(Fcompvec));
end

%Fvec = zeros(nvox,Nk);
Fvec = zeros(nvox, prod(Nk));


for i=1:nvox           
    
    %Fvec(i,:) = sum(repmat(weights(i,:),[Nk 1]).*Fcompvec, 2);    
    Fvec(i,:) = sum(repmat(weights(i,:),[prod(Nk) 1]).*Fcompvec, 2);
    
    %Fvecloop(i,:) = sum(bsxfun(@times,weights(i,:), Fcompvec), 2);    
end


%tic;
%Fvecloop = sparse(Nk1,Nk2);
%for i=1:ncomp    
%     Fvecloop = Fvecloop + weights(:,i).*Fcompvecarray(:,i);
%end    
%toc;


%construct matrix where every row is K*Fvec(z_n) for voxel n

%KalphaFvecloop = zeros(nvox,nvoxaug);
% tic;
% for i=1:nvox
%     KalphaFvecloop(i,:) = Kalpha*Fvec(i,:)';
% end
% toc;

%KalphaFvec = (Kalpha * Fvec')';    
KalphaFvec = (sparse(Kalpha) * sparse(Fvec)')';   



%calculate the log-likelihood for all voxels
if options.ILT.reg   
    
    LogLi = (-1./(2*sigvec.^2)) .* sum((allimgaug - KalphaFvec).^2,2);

else            
    
    LogLi = (-1./(2*sigvec.^2)) .* sum((allimg - KalphaFvec).^2,2);
     
     
    RESNORM = norm(allimg-KalphaFvec);
end



end