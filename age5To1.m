function [valA1] = age5To1(valA5)
[ageDim, valDim] = size(valA5);
valA1 = zeros(ageDim*5 , valDim);

for i = 1 : ageDim
    valA1(((i-1)*5+1) : i*5 , :) = ones(5 , valDim) .* valA5(i , :);
end
    




         
         