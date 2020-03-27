function S = KernelDK(params,gradechoinv)
%kernel for isotropic diffusivity with kurtosis

[~, b, ~, ~, ~] = unpack_gradechoinv(gradechoinv);

d = params(1);
k = params(2);

S = exp (-b.*d + b.^2*d^2*k/6 );

