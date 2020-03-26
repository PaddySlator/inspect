function Vf = integrate_spectrum(S,T1,T2,T1bounds,T2bounds)

%T1bounds and T2bounds are 

if isempty(T2bounds)
    Vf = zeros(1, length(T1bounds) - 1);
    for i=1:(length(T1bounds) - 1)
        
        Vf(i) = sum(sum( S(T1bounds(i) < T1 & T1 < T1bounds(i+1),:) ));
        
    end
    Vf = Vf/sum(Vf);
end


if isempty(T1bounds)
    Vf = zeros(1, length(T2bounds) - 1);
    for i=1:(length(T2bounds) - 1)
        
        Vf(i) = sum(sum( S(T2bounds(i) < T2 & T2 < T2bounds(i+1),:) ));
        
    end
    Vf = Vf/sum(Vf);
end




% figure;
% contour(T2,T1,S);
% set(gca,'yscale','log')







end
