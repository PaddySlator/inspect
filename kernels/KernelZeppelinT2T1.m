function S = KernelZeppelinT2T1(params,protocol)

dpar = params(1);
dperp = params(2);
theta = params(3);
phi = params(4);
t2 = params(5);
t1 = params(6);


S = Zeppelin([dpar dperp theta phi],protocol) .* T2decay(t2,protocol) .* T1inv(t1,protocol);
