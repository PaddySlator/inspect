function S = KernelMeas(kernel,gradechoinv)

%and returns a function handle for that kernel

%protocol gives the experimental  parameters  which  are  varied  to  yield  contrast in intrinsic MR properties w.





fun = str2func(['Kernel' kernel.name]);

S = fun(kernel.params,gradechoinv);








        

end
        








