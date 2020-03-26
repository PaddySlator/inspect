function S = KernelDT2T1(params,gradechoinv)

d = params(1);
t2 = params(2);
t1 = params(3);



S = KernelD(d,gradechoinv) .* KernelT2(t2,gradechoinv) .* KernelT1inv(t1,gradechoinv);




