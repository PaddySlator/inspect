function S = KernelT2(params,gradechoinv)

%kernel for t2 decay


[~, ~, TE, ~, ~] = unpack_gradechoinv(gradechoinv);

t2 = params(1);



S = exp( -TE ./ t2);
       




