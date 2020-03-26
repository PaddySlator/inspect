function S = KernelT1inv(params,gradechoinv)

%kernel for t1 inversion recovery

%calculate the signal for given contrast encodings and MR contrast value

%t are experimental  parameters  which  are  varied  to  yield  contrast in intrinsic MR properties w.

% t - [TI TE TR]
% w - T1 value 

[~, ~, ~, TI, TR] = unpack_gradechoinv(gradechoinv);

t1 = params(1);


%calculate the signal                    
S = abs(1 - 2 * exp( - TI / t1 ) + exp(-TR/t1) );

                    
