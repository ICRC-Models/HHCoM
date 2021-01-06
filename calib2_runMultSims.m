% Parameter estimation

% Calculates the negative log likelihood of a sample of potential parameter sets
% Accepts:
% 1) paramSetIdx (index in larger sample matrix to start selection of sub-set)
% Saves:
% 1) File: negSumLogL_calib_[date].dat (negative log likelihood for each parameter set in sub-set)

function calib2_runMultSims(paramSetIdx , tstep_abc , date_abc)

t_curr = tstep_abc;
date = date_abc;

%% Cluster information
pc = parcluster('local');    % create a local cluster object
pc.JobStorageLocation = strcat('/gscratch/csde/carajb' , '/' , getenv('SLURM_JOB_ID'))    % explicitly set the JobStorageLocation to the temp directory that was created in the sbatch script
numCPUperNode = str2num(getenv('SLURM_CPUS_ON_NODE'))
parpool(pc , numCPUperNode)    % start the pool with max number workers

nPrlSets = 5; %numCPUperNode;

%% Load all particles
paramDir = [pwd , '/Params/'];
masterSetMatrix = load([paramDir , 'masterSets_calib_' , date , '_' , num2str(t_curr) , '.dat']); % load most recent parameter sample
pIdx = load([paramDir , 'pIdx_calib_' , date , '_0' , '.dat']); % load parameter indices
orderedLL = load([paramDir , 'orderedLL_calib_' , date , '_' , num2str(t_curr) , '.dat']); % load most recent ordered log-likelihoods

%% Get indices and parameter values of top 25 sets (*(25/10220))
numBestFits = 25;
numSets25 = size(orderedLL , 1)*(numBestFits/10220);
top25Inds = orderedLL(1:numSets25,1);
top25Params = masterSetMatrix(:,top25Inds);

%% If on Phase 2 of calibration, uncomment the following to plot resampled Ph1 parameters + Ph2 parameters
pIdx = load([paramDir,'pIdx_calib_' , date , '_0_wPh1Resample.dat']);
masterResampleSubsetMatrix = load([paramDir , 'masterResampleSubsetMatrix_calib_' , date , '_' , num2str(t_curr) , '.dat']); % load most recent Ph1 resampled parameters
masterCombinedPhaseMatrix = [masterResampleSubsetMatrix ; masterSetMatrix];
top25Params = masterCombinedPhaseMatrix(:,top25Inds);

%% Set up paramsSub for indexing into paramSet matrix
[paramsAll] = genParamStruct();
paramsSub = cell(length(pIdx),1);
startIdx = 1;
for s = 1 : length(pIdx)
    paramsSub{s}.length = paramsAll{pIdx(s)}.length;
    paramsSub{s}.inds = (startIdx : (startIdx + paramsSub{s}.length - 1));
    startIdx = startIdx + paramsSub{s}.length;
end

%% Obtain model output for each set of sampled parameters
%for m = 1 : nPrlSets : numBestFits
    subMatrixInds = [paramSetIdx : (paramSetIdx + nPrlSets - 1)];
    %subMatrixInds = [m : (m + nPrlSets - 1)];
    parfor n = 1 : nPrlSets
        paramSet = top25Params(:,subMatrixInds(n));
        futureSim(1 , pIdx , paramsSub , paramSet , (paramSetIdx + n - 1) , tstep_abc , date_abc);
        %historicalSim(1 , pIdx , paramsSub , paramSet , (paramSetIdx + n - 1) , tstep_abc , date_abc);
    end
%end

