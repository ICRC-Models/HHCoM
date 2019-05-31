% Identify x-many best-fitting parameter sets

numSets = 25;

paramDir = [pwd , '/Params/'];
negSumLogLmatrix = load([paramDir , 'negSumLogL_calib_29May19.dat']);
paramSetMatrix = load([paramDir , 'paramSets_calib_29May19.dat']);

numSubsets = size(negSumLogLmatrix,1)/17;

negS_format = reshape(negSumLogLmatrix , [17,numSubsets]);
maxV = max(negS_format(1,:));
setVec = [1:16:maxV];
for j = 1 : length(setVec)
     if ~any(setVec(j) == negS_format(1,:)) 
         negS_format(:,end+1) = [setVec(j); ones(16,1).*10000000]; 
     end
end

[temp,firstRowOrder] = sort(negS_format(1,:));
negS_ordered = negS_format(:,firstRowOrder);

negS_ordered_flatDat = reshape(negS_ordered(2:end,:),[16*numSubsets,1]);
[vals,inds] = sort(negS_ordered_flatDat,'ascend');

file = 'bestParamSets_calib_29May19.dat';
paramDir = [pwd , '/Params/'];

for i = 1:numSets
    dlmwrite([paramDir , file] , [vals(i); inds(i); paramSetMatrix(:,i)] , 'delimiter' , ',' , ...
    'roffset' , 1 , 'coffset' , 0 , '-append')
end
