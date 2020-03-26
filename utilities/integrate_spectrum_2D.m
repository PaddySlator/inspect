function Vf = integrate_spectrum_2D(output,w1bounds,w2bounds,ploton)

%T1bounds, T2bounds and T3bounds are matrices where each row defines a
%region to integrate the spectrum within, 
% so if the first rows are: T1bounds = [0 10],T2bounds = [0 20],T3bounds = [50 100],
% then the first region to integrate in is 0 < T1 < 10, 0 < T2 < 20, 50 <
% T3 < 100


F = output.F;
w1 = output.w1;
w2 = output.w2;
%w3 = output.w3;

%F = MERA_fit.MERA_out2D.S;
%w1 = MERA_fit.MERA_out2D.T;
%w2 = MERA_fit.MERA_out2D.T2;


%number of compartments to integrate over
% Ncomp1= length(T1bounds) - 1;
% 
% Vf = zeros(Ncomp1,1);
% 
% for i=1:Ncomp1    
%     Vf(i) = sum(sum(sum(F(T1bounds(i) < T1 & T1 < T1bounds(i+1), :, :) )));
% end
% 
% Vf = Vf/sum(Vf);
% 
% 
% Ncomp2 = length(T2bounds) - 1;
% for i=1:Ncomp2
%     %Vf(i) = sum(sum(sum(F(:, T2bounds(i) < T2 & T2 < T2bounds(i+1), :) )));    
% end


if ~(size(w1bounds,2)==2 && size(w1bounds,2)==2 && size(w1bounds,2)==2)
    error('need 2 column matrices to define the upper and lower limits of the boundaries')
end

Nw1 = size(w1bounds,1);
Nw2 = size(w2bounds,1);
%Nw3 = size(w3bounds,1);

% if range([Nw1 Nw2 Nw3]) ~= 0
%     error('the boundary matrices should be the same size')   
% end


for i=1:Nw1
    Vf(i) = sum(sum(F(w1bounds(i,1) <= w1 & w1 <= w1bounds(i,2),...
        w2bounds(i,1) <= w2 & w2 <= w2bounds(i,2))));
end


% for i=1:Nw1    
%     Vf(i) = sum(sum(sum(F(w1bounds(i,1) <= w1 & w1 <= w1bounds(i,2),...
%         w2bounds(i,1) <= w2 & w2 <= w2bounds(i,2),...
%         w3bounds(i,1) <= w3 & w3 <= w3bounds(i,2)))));
% end
Vf = Vf/sum(Vf);










    

    


% figure;
% contour(T2,T1,S);
% set(gca,'yscale','log')







end
