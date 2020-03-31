function inspect_sim = simulate_inspect_map(spectra,vfmap,gradechoinv,kernel)

%Forward Model for inspect continuous mapping version.
%Given canonicial basis spectra, and ground truth volume fraction maps
%assumes that the signal from each voxel is a weighted sum of 
%contributions from a set of canonicial basis spectra

% REQUIRED INPUTS:
%
% spectra - canonical spectral components 
%
% vfmap - image containing the ground truth volume fractions 
%
% gradechoinv - array, or path to a .txt file, of the MR acquisition 
% parameters.  Format: [gx gy gz b TE TI TR] where gx gy gz are the gradient 
% directions, b is the b-value, TE is the echo time, TI is the inversion 
% time, and TR is the repetition time. 
%
% kernel - string specifying the choice of kernel, which relates
% the MR acquisition parameters to the MR signal. For example:
% - 'DT2' for diffusion-T2, for which the kernel equation 
%   is exp(-b*D)*exp(-TE/t2). 
% - 'T1' for T1 inversion recovery, for which the kernel equation 
%   is abs(1-2*exp(-TI/t1) + exp(-TR/t1)).
%
% OPTIONAL INPUTS:
%
% options - algorithm options
%
%
%
%
% Author: Paddy Slator, p.slator@ucl.ac.uk
%
%
% LICENSE
% <inspect toolbox for qMRI analysis>
% Copyright (C) <2020>  <Paddy J. Slator, p.slator@ucl.ac.uk>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.







