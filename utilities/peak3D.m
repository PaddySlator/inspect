function [Av,Fv,Tv,Nc] = peaks3D(S)
% automatic 3D peak identification

[r,c] = size(S);
zz = find(S); xz = min(S(zz));

%          zz = S; xz = min(S(zz));

[Si,Nc] = bwlabel(im2bw(S./xz));


Tv = zeros(Nc,2);
Av = zeros(Nc,1);

for kc = 1:Nc
    [rc,cc] = find(Si == kc);
    fc = diag(S(rc,cc));
    Tv(kc,1) = exp(sum(fc.*log(fitting.T(rc)))./sum(fc)); %geometric mean
    Tv(kc,2) = exp(sum(fc.*log(fitting.T2(cc)))./sum(fc)); %geometric mean
    Av(kc) = sum(fc);
end
%volume fraction of each peak
Fv=Av/sum(Av);

            
            



end