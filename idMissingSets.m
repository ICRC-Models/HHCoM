function idMissingSets(tstep_abc , date_abc , nSets)

t_curr = tstep_abc;
date = date_abc;

%% Load all parameter sets
paramDir = [pwd , '/Params/'];
paramSetMatrix = load([paramDir , 'paramSets_calib_' , date , '_' , num2str(t_curr) , '.dat']); % load most recent parameter sample
negSumLogLmatrix = load([paramDir , 'negSumLogL_calib_' , date , '_' , num2str(t_curr) , '.dat']); % load most recent log-likelihoods

%% Filter out failed parameter sets (timed-out, etc.)
numSubsets = size(negSumLogLmatrix,1)/5; % calculate number of sub-sets that actually ran (vs. timed-out, failed, etc.)
negS_format = reshape(negSumLogLmatrix , [5,numSubsets]); % first row is paramSetIdx, next 4 rows log-likelihoods for that sub-set
setVec = [1:4:nSets];
missingV = [];
for j = 1 : length(setVec) % identify failed parameter sets
     if ~any(setVec(j) == negS_format(1,:))
         missingV = [missingV , setVec(j)];
     end
end

fileF = ['missingSets_calib_', date , '_' , num2str(t_curr) , '.dat'];
dlmwrite([paramDir, fileF] , missingV, 'delimiter',' ');

