function scatterpoints = spectrum_to_scatter(output)
%go from spectrum to points for a scatter plot

F = output.F;
w1 = output.w1;
w2 = output.w2;
w3 = output.w3;

Nw1 = length(output.w1);
Nw2 = length(output.w2);
Nw3 = length(output.w3);


scatterpoints = [];

for i=1:Nw1
    for j=1:Nw2
        for k=1:Nw3
            if F(i,j,k)~=0 
                scatterpoints = [scatterpoints; w1(i) w2(j) w3(k)];
            end
        end
    end
end



end