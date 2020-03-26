function F = construct_spectrum_from_components(Fcomp,weights)

% Fcomp are the spectral components 
% z is a M by 1 vector with the weights of the component to the spectrum

F = zeros(size(Fcomp{1}));

for i=1:length(Fcomp)  
    F = F + weights(i)*Fcomp{i};  
end



%normalise?



end




