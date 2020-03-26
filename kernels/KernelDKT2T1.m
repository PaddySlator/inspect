function S = KernelDKT2T1(params,gradechoinv)

d = params(1);
k = params(2);
t2 = params(3);
t1 = params(4);



S = KernelDK([d k],gradechoinv) .* KernelT2(t2,gradechoinv) .* KernelT1inv(t1,gradechoinv);




