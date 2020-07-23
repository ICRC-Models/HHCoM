function [paramSet] = even_sample(boundsList,valuesPerBound)
%even_sample creates an evenly spaced multidimensional sample
%   boundslist is a (#bounds by 2) matrix. Each row has a lower and upper
%   bound for a given parameter.
%   valuesPerBound is the sampling density (multidimensional grid cell
%   length)
%   returns the multidimensinal sample as a #bounds by
%   valuesPerBound^#bounds matrix
numBounds = size(boundsList,1);
numTrials = valuesPerBound.^numBounds;
mat = zeros(numBounds,numTrials);
counters = zeros(numBounds,1);
rangesizes = boundsList(:,2)-boundsList(:,1);

workingIdx = 1;
recipVPB = 1/(valuesPerBound-1);%3 samples calls for low bound, midpoint, and highpoint
matCellCounter = 1;

while matCellCounter <= numTrials;
    
    disp(counters);
    mat(:,matCellCounter) = boundsList(:,1) + counters.*rangesizes.*recipVPB;
    matCellCounter = matCellCounter+1;
    
    if counters(1)<valuesPerBound-1
        counters(1)=counters(1)+1;
    else
        searching = 1;
        while(searching)
            
            workingIdx = workingIdx+1;
            if workingIdx>numBounds
                searching =0;
            elseif counters(workingIdx)<valuesPerBound-1
                counters(workingIdx) = counters(workingIdx)+1;
                counters(1:workingIdx-1) = zeros(workingIdx-1,1);
                workingIdx =1;
                searching =0;
            end
        end
           
end


paramSet =  mat;
end

