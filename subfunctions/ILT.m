function output = ILT(sig,gradechoinv,options)

%Dicretized inverse Laplace transform for arbitrary kernal and contrast
%encondings 

%solve S = KF for F
% S is the signal  
% K is the matrix of kernals
% F is the spectrum  


%gradechoinv - values of the MR contrast encodings
%e.g. gradient direction, b-values, echo times, inversion time



% options - algorithm options
% including options.kernel - 
%possible kernels



%unpack the contrast encodings
%[g, b, TE] = unpack_gradechoinv(gradechoinv);


%make camino protocol
%protocol = gradechoinv2protocol(gradechoinv);



%get the number of measurements 
Nmeas = size(gradechoinv,1);
%Nmeas = protocol.totalmeas;

%define the grid on which to estimate the diffusion/relaxometry spectrum
Nk = options.Nk;
mink = options.mink;
maxk = options.maxk;


%Nk1 = options.Nk1;
%Nk2 = options.Nk2;
%Nk3 = options.Nk3;
% 
% for i=1:length(Nk)
%     if options.loggrid(i)
%         w(:,i) = logspace(log10(mink(i)), log10(maxk(i)),Nk(i));  
%     else
%         w(:,i) = linspace(mink(i), maxk(i), Nk(i));
%     end
% end

%w1 = logspace(log10(options.mink1),log10(options.maxk1),Nk1);
%w2 = linspace(options.mink2,options.maxk2, Nk2);
%w3 = logspace(log10(options.mink3),log10(options.maxk3),Nk3);


%build the signal matrix - Nmeas x 1 
if iscolumn(sig)
    S = sig;
else
    S = sig';
end




%build kernel function - takes contrast encodings and discretisation of the grid
%as inputs


%if given as input
if isfield(options,'Kalpha')
    Kalpha = options.Kalpha;
    K = options.K;
    grid = options.grid;
elseif isfield(options,'K')
    K = options.K;
    grid = options.grid;
else
    %build the matrix of kernals - Nmeas x prod(Nk)
    %K = zeros(Nmeas,Nk1*Nk2*Nk3);
    K = zeros(Nmeas, prod(Nk));
    
    grid = getkernelgrid(options);
    
    allgridcombs = ( combvec(grid{:}) )';
    Ngridcombs = size(allgridcombs,1);
    
    kernel.name = options.kernel;
    
    for i=1:Ngridcombs
        kernel.params = allgridcombs(i,:);
        K(:,i) = KernelMeas(kernel,gradechoinv);
    end
            
    if options.reg
        %augment kernal matrix
        H = speye(Ngridcombs, Ngridcombs);
        %H(end,end) = 0;
        Kalpha = [K; options.alpha*H];      
    end
end


%option to not actually fit the ILT - just get the dictionary values
if isfield(options,'onILT')
    if ~options.onILT 
        output.grid = grid;
        output.S = S;
        output.K = K;
        if options.reg
            output.Kalpha = Kalpha;
        end
        return
    end
end
    
% Kloop = zeros(Nmeas,Nk1*Nk2*Nk3);


% for meas = 1:Nmeas
%     l=1;
%     for k=1:Nk3
%         for j=1:Nk2
%             for i=1:Nk1         
%                 Kloop(meas,l) = exp(-b(meas)*w1(i)) * ...
%                     exp(-TE(meas)/w2(j)) * ...
%                     abs(1 - 2*exp(-(TI(meas) + TE(meas))/w3(k)) + exp(-TR(meas)/w3(k)) );
%                     ...abs(1 - 2*exp(-(b3(meas))/w3(k)) + exp(-TR/w3(k)) );
%                     %abs(1 - 2*exp(-b3(meas)/w3(k)) );
%                 l=l+1;
%             end
%         end
%     end
% end
% output.Kloop = Kloop;







if options.reg     
    %augment data vector
    Ngridcombs = size(K,2);          
    H = speye(Ngridcombs, Ngridcombs);
    S = [S; zeros(size(H,1),1)];
    %nnls with regularisation   
    [F,RESNORM] = lsqnonneg(Kalpha,S);    
    %possible faster way
    %F = tntnn(K,S,alpha);       
    %RESNORM = norm(S - K*F);
    %another possible faster way
    %[F,RESNORM] = fnnls(Kalpha'*Kalpha, Kalpha'*S);
    %and another!
    %[F,RESNORM] = nnls(Kalpha,S);
else
    %nnls without regularisation
    [F,RESNORM] = lsqnonneg(K,S);  
    %possible faster way    
    %F = tntnn(K,S);       
    %RESNORM = norm(S - K*F);
    %another possible faster way
    %[F,RESNORM] = fnnls(K'*K, K'*S);
    %and another!
    %[F,RESNORM] = nnls(K,S);
end


%store the vectorised spectrum
Fvec = F;
%reshape F if 2D or above
if ~isscalar(Nk)
    F = reshape(F, Nk);
end


output.grid = grid;

output.S = S;
output.F = F;
output.Fvec = Fvec;
output.RESNORM = RESNORM;
%output.allgridcombs = allgridcombs;
output.gradechoinv = gradechoinv;
output.K = K;
if ~options.reg %equivalent to regularised case with alpha=0 - return so that downstream functions can just use Kalpha as input
    Kalpha = K; 
end
output.Kalpha = Kalpha;






