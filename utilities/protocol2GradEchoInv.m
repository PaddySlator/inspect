function GradEchoInv = protocol2GradEchoInv(protocol,saveon)
% Takes camino.m protocol file and converts to "gradechoinv file"
%
% gradechoinv is a file containing diffusion-weightings, echo times and inversion
% times for a combined diffusion-relaxometry acquisition. Each row is [x y z b TE TI], 
%   - x, y, z are the diffusion gradient directions
%   - b is the b-value
%   - TE is the echo time 
%   - TI is inversion time 
% The size of gradechoinv is the total number of measurements times 6


GradEchoInv = nan(protocol.totalmeas,6);

%gradient directions
GradEchoInv(:,1:3) = protocol.grad_dirs;

%get the b-values
GradEchoInv(:,4) = GetBvalues_smm2(protocol);

%get the echo times
if isfield(protocol,'TE')
    GradEchoInv(:,5) = protocol.TE;
end

%get the inversion times
if isfield(protocol,'TI')
    GradEchoInv(:,6) = protocol.TI;
end

if saveon
   fid = fopen('grad_echo_inv.txt','w');
   fprintf(fid, '%f %f %f %f %f %f \n', GradEchoInv')
end


end