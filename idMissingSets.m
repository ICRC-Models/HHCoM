%% Load all parameter sets
paramDir = [pwd , '/Params/'];
paramSetMatrix = load([paramDir , 'paramSets_calib_' , '18July19_4', '.dat']); % load most recent parameter sample
negSumLogLmatrix = load([paramDir , 'negSumLogL_calib_' , '18July19_4', '.dat']); % load most recent log-likelihoods

%% Filter out failed parameter sets (timed-out, etc.)
numSubsets = size(negSumLogLmatrix,1)/17; % calculate number of sub-sets that actually ran (vs. timed-out, failed, etc.)
negS_format = reshape(negSumLogLmatrix , [17,numSubsets]); % first row is paramSetIdx, next 16 rows log-likelihoods for that sub-set
setVec = [1:16:3114];
missingV = [];
for j = 1 : length(setVec) % identify failed parameter sets
     if ~any(setVec(j) == negS_format(1,:))
         missingV = [missingV , setVec(j)];
     end
end

fileF = ['missingSets_calib_', '18July19_4', '.dat'];
dlmwrite([paramDir, fileF] , missingV, 'delimiter',' ');

