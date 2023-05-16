function S = KernelT2T1(params,gradechoinv)

t2 = params(1);
t1 = params(2);

S = KernelT2(t2,gradechoinv) .* KernelT1inv(t1,gradechoinv);




