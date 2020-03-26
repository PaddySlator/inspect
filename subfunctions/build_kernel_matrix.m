function K = build_kernel_matrix(kernels,gradechoinv,options)

%options for kernels are 
% expdecay   - exponential decay, i.e. exp(-tw)
% expdecayrecip  - exponential decay with reciprocal parameter, i.e. exp(-t/w)
% invrecov - invesion recovery, exp(1-...)

Nkernals = length(kernels);

expdecay = strcmp(kernels,'expdecay');
expdecayrecip = strcmp(kernels,'expdecayrecip');
invrecov = strcmp(kernels,'invrecov');



Nk = options.Nk;

loggrid = [1 0];

Ndim = length(Nk);

Nmeas = size(gradechoinv,1);

%build the grid for calculating the kernal on
for i=1:Ndim
    if loggrid(i)
        w(:,i) = logspace(log10(options.mink(i)),log10(options.maxk(i)),Nk(i));        
    else
        w(:,i) = linspace(options.mink(i),options.maxk(i), Nk(i));
    end
end

%get the sampling times 
%b-value
t1 = gradechoinv(:,4);
%echo time
t2 = gradechoinv(:,5);

t = gradechoinv(:,4:5);

%build the matrix of kernals - i.e. the dictionary - size Nmeas x w1*w2*w3

K = zeros(Nmeas,prod(Nk));

for meas = 1:Nmeas
    l=1;
    for j=1:Nk(2)
        for i=1:Nk(1)
          
            
            K(meas,l) = sum(expdecay.*[exp(-t1(meas)*w(i,1))    exp(-t2(meas)*w(i,2))]) * ...
                sum(expdecayrecip.*[exp(-t1(meas)/w(j,1))   exp(-t2(meas)/w(j,2))]);
            l=l+1;
        end
    end
end











end