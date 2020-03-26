function [c, ceq] = probcon(p)
%function which constrains parameter vector

%bounding nonlinear constraint - not required
c = [];
%constrain this to be zero
ceq = sum(p) - 1;


end



