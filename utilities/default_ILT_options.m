function options = default_ILT_options(kernelname,gradechoinv)
%default options for the ILT given the kernel and acquisition protocol

options.kernel = kernelname;


%default no regularisation - faster
%options.reg = 1;
%options.alpha = 0.01;
options.reg = 0;
options.alpha = 0;




%get the parameter names for this kernel       
paramstrings = GetKernelParameterStrings(kernelname);
nparam = length(paramstrings);


%size of the grid on which to estimate the spectrum 
if nparam == 1
    options.Nk = 200;
elseif nparam == 2
    options.Nk = [50 50];
elseif nparam == 3
    options.Nk = [20 20 20];    
elseif nparam == 4
    options.Nk = 10 * ones(1,nparam);    
else
    options.Nk = 10 * ones(1,nparam);  
end
    

%minimum and maximum values for the grid for each parameter
%diffusivity

%get reasonable max and min values from the acquisition protocol

[g, b, TE, TI, TR] = unpack_gradechoinv(gradechoinv);

%get the unique non-zero b-values
b = unique(b(b~=0));
%heuristic choice of how far above/below the 1/min(b) to have the max/min D
bcoeff = 5;
min_grid.d = (1/bcoeff) * 1/max(b);
max_grid.d = bcoeff * 1/min(b);

min_grid.dpar = (1/bcoeff) * 1/max(b);
max_grid.dpar = bcoeff * 1/min(b);

min_grid.dperp = (1/bcoeff) * 1/max(b);
max_grid.dperp = bcoeff * 1/min(b);

min_grid.d1 = (1/bcoeff) * 1/max(b);
max_grid.d1 = bcoeff * 1/min(b);

min_grid.d2 = (1/bcoeff) * 1/max(b);
max_grid.d2 = bcoeff * 1/min(b);



%min_grid.d = 2*10^-4;
%max_grid.d = 5;
%max_grid.d = 0.005;
%this needs to be lower if we are including kurtosis - otherwise dictionary
%of signals goes to infinity - could improve by having "non-square"
%dictionary
if any(strcmpi(paramstrings,'k'))
    max_grid.d = 0.003;
end

%kurtosis 
min_grid.k=0;
max_grid.k=3;

%T2
TE = unique(TE);
%heuristic choice of how far above/below the min/max TE to have the max/min
%T2
TEmincoeff = 5;
TEmaxcoeff = 5;

min_grid.t2 = (1/TEmincoeff) * min(TE);
max_grid.t2 = TEmaxcoeff * max(TE);

%min_grid.t2 = 5;
%max_grid.t2 = 200;

%min_grid.t2 = 0.005;
%max_grid.t2 = 0.2;


%T1
%heuristic choice of how far above/below the min/max TI to have the max/min
%T1
TI = unique(TI);

TImincoeff = 1;
TImaxcoeff = 1;

min_grid.t1 = (1/TImincoeff) * min(TI);
max_grid.t1 = TImaxcoeff * max(TI);

%min_grid.t1 = 250;
%max_grid.t1 = 5000;

%min_grid.t1 = 0.25;
%max_grid.t1 = 5;


%fill the min/max grid options for each parameter
for i = 1:nparam
    options.mink(i) = min_grid.(paramstrings{i});
    options.maxk(i) = max_grid.(paramstrings{i});          
end

% choose whether the grid points are log-scaled. Choose log scale if peaks
% are expected to have large (i.e. order of magnitude) separatation. 

%Log-scale for d - e.g. diffusivity and pseudo-diffusivity
logscale.d = 1;
logscale.dpar = 1;
logscale.dperp = 1;
logscale.d1 = 1;
logscale.d2 = 1;
%linear scale for t1/t2 as these are typically not widely separated (at least
%for the experiments/protocols thus far)
logscale.t2 = 0;
logscale.t1 = 0;

logscale.k = 0;


%fill the log-scale options for each parameter
for i = 1:nparam
    options.loggrid(i) = logscale.(paramstrings{i});
end

