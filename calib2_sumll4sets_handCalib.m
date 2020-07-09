% Parameter estimation

% Run simulations with trial parameter sets in parallel and calculate the
% corresponding negative summed-log-likelihoods for model outcomes compared to observed data.
% Accepts:
% 1) paramSetIdx (column index into larger parameter trial matrix
% (paramSetMatrix) to start selection of sub-set). Should equal 1 unless the
% number of columns in paramSetMatrix > 12, or the maximum number of trial
% sets you can run in parallel on the Simulation cluster
% Saves:
% 1) File: negSumLogL_calib_[date]_[t_curr].dat (negative log likelihood for each parameter set in sub-set)

function calib2_sumll4sets_handCalib(paramSetIdx , tstep_abc , date_abc)

t_curr = tstep_abc; % calibration iteration
date = date_abc; % date

createParamSet_handCalib(tstep_abc , date_abc); % create files with specified pIdx and paramSetMatrix values

%% Load parameters
paramDir = [pwd ,'/Params/'];
paramSetMatrix = load([paramDir,'paramSets_calib_' , date , '_' , num2str(t_curr) , '.dat']); % load matrix of trial parameter values (size: [# parameters x # trials])
nPrlSets = 3; % # of trials, or different parameter sets you want to run in parallel (should be <= 12) 
subMatrixInds = [paramSetIdx : (paramSetIdx + nPrlSets - 1)]; % list of indices for the columns in paramSetMatrix holding the parameter sets you will run in parallel

pIdx = load([paramDir,'pIdx_calib_' , date , '_' , num2str(tstep_abc) , '.dat']);  % load indices into paramsAll cell array for parameters you want to include in calibration

%% Set up paramsSub for indexing into paramSet matrix
[paramsAll] = genParamStruct(); % load cell array holding information about all parameters that are able to be calibrated
paramsSub = cell(length(pIdx),1); % create an empty cell array to hold relevant information about subset of parameters included in this calibration
startIdx = 1;
for s = 1 : length(pIdx) % loop through parameters included in calibration
    paramsSub{s}.length = paramsAll{pIdx(s)}.length;
    paramsSub{s}.inds = (startIdx : (startIdx + paramsSub{s}.length - 1)); % save indices into the rows of paramSetMatrix for the relevant parameters
    startIdx = startIdx + paramsSub{s}.length;
end

%% Obtain model output for each trial parameter set
negSumLogLSet = zeros(nPrlSets,1);
parfor n = 1 : nPrlSets % start nPrlSets-many historical model simulations in parallel
    paramSet = [paramSetMatrix(:,subMatrixInds(n))]; % column in paramSetMatrix with the trial parameters to use for a simulation
    %futureSim(1 , pIdx , paramsSub , paramSet , (paramSetIdx + n - 1) , tstep_abc , date_abc);
    [negSumLogL] = historicalSim(1 , pIdx , paramsSub , paramSet , (paramSetIdx + n - 1) , tstep_abc , date_abc); % start historical simulation using the trial parameter set (paramSet)
    negSumLogLSet(n,1) = negSumLogL;
end

%% Save parameter set indices and their corresponding negSumLogL values
formatOutput = zeros(nPrlSets+(nPrlSets/4) , 1);
paramSetIdxBY4 = [paramSetIdx : 4 : (paramSetIdx + nPrlSets - 1)];
nPrlSetsBY4 = [1 : 4 : (nPrlSets - 1)];
for i = 1 : (nPrlSets/4)
    formatOutput(nPrlSetsBY4(i) + (i-1) : nPrlSetsBY4(i) + 4 + (i-1) , 1) = [paramSetIdxBY4(i) ; negSumLogLSet(((i-1)*4+1 : i*4) , 1)];
end

file = ['negSumLogL_calib_' , date , '_' , num2str(t_curr) , '.dat'];
paramDir = [pwd , '/Params/'];
dlmwrite([paramDir, file] , formatOutput , 'delimiter' , ',' , 'roffset' , 1 , 'coffset' , 0 , '-append' , 'precision' , 9)

% for j = 1 : nPrlSets
%     pathModifier = ['toNow_' , date , '_noBaseVax_baseScreen_hpvHIVcalib_' , num2str(t_curr) , '_' , num2str(paramSetIdx + j - 1)];
%     savdir = [pwd , '/HHCoM_Results/'];
%     if negSumLogLSet(j,1) < -2000 %-200000.00
%         delete([savdir , pathModifier , '.mat']);
%     end
% end
