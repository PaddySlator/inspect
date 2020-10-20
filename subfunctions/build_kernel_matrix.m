function K = build_kernel_matrix(gradechoinv,options)

% 

% Nkernals = length(kernels);
% 
% expdecay = strcmp(kernels,'expdecay');
% expdecayrecip = strcmp(kernels,'expdecayrecip');
% invrecov = strcmp(kernels,'invrecov');



Nk = options.Nk;

Nmeas = size(gradechoinv,1);


%build the matrix of kernals - Nmeas x prod(Nk)
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