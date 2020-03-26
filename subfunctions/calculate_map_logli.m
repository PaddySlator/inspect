function LogLi = calculate_map_logli(S,weights,Fcomp,sig,Kalpha)

F = construct_spectrum_from_components(Fcomp,weights);

Fvec = F(:);



%LogLi = sum( - (1/(2*sig^2))*(S - Kalpha*Fvec).^2);


LogLi =  - (1/(2*sig^2)) * sum( (S - Kalpha*Fvec).^2 );

end