function grid = getkernelgrid(options)

%define the grid on which to estimate the diffusion/relaxometry spectrum
%This depends on the choice of kernel, and the number of gridpoints in each
%dimension

paramstrings = GetKernelParameterStrings(options.kernel);
Nparam = length(paramstrings);

Nk = options.Nk;
mink = options.mink;
maxk = options.maxk;


for i=1:length(Nk)
    if options.loggrid(i)
        grid{i} = logspace(log10(mink(i)), log10(maxk(i)),Nk(i));  
    else
        grid{i} = linspace(mink(i), maxk(i), Nk(i));
    end
end




% k1min = options.mink1; 
% k1max = options.maxk1;
% k2min = options.mink2;
% k2max = options.maxk2;
% k3min = options.mink3;
% k3max = options.maxk3;
% 
% anglemin = 0;
% anglemax = 2*pi;
% 
% 
% Nk1 = options.Nk1;
% Nk2 = options.Nk2;
% Nk3 = options.Nk3;
% 
% switch options.kernel
%     case 'BallT2'
%     
%     case 'BallT2T1'
%         %grid(:,1) = logspace(log10( dmin ),log10( dmax ), Nk1);
%         %grid(:,2) = linspace( t2min , t2max , Nk2);
%         %grid(:,3) = logspace(log10( t1min ),log10( t1max ),Nk3);
%         
%         grid{1} = logspace(log10( k1min ),log10( k1max ), Nk1);
%         grid{2} = linspace( k2min , k2max , Nk2);
%         grid{3} = logspace(log10( k3min ),log10( k3max ),Nk3);
%         
%     case 'ZeppelinT2T1'
%         grid(:,1) = logspace(log10( dmin ),log10( dmax ), Nk1);
%         grid(:,2) = logspace(log10( dmin ),log10( dmax ), Nk1);
%         grid(:,3) = linspace(anglemin, anglemax, Nk1);
%         grid(:,4) = linspace(anglemin, anglemax, Nk1);
%         grid(:,5) = linspace( t2min , t2max , Nk1);
%         grid(:,6) = logspace(log10( t1min ),log10( t1max ),Nk1);
%                       
%         
%         
%         
% end