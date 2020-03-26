function S = KernelDT2(params,gradechoinv)

d = params(1);
t2 = params(2);

S = KernelD(d,gradechoinv) .* KernelT2(t2,gradechoinv);



