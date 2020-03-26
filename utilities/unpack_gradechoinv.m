function [g, b, TE, TI, TR] = unpack_gradechoinv(gradechoinv)

%unpack the contrast encodings
g = gradechoinv(:,1:3);
b = gradechoinv(:,4);
if size(gradechoinv,2) >= 5
    TE = gradechoinv(:,5);
else
    TE = NaN;
end
if size(gradechoinv,2) >= 6
    TI = gradechoinv(:,6);
else
    TI = NaN;
end
if size(gradechoinv,2) >= 7
    TR = gradechoinv(:,7);
else
    TR = NaN;
end




end