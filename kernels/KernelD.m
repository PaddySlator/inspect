function S = KernelD(params,gradechoinv)
%kernel for isotropic diffusivity

[~, b, ~, ~, ~] = unpack_gradechoinv(gradechoinv);

d = params(1);

S = exp (-b .* d);




