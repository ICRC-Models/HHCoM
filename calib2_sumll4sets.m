% Parameter estimation

% Calculates the negative log likelihood of a sample of potential parameter sets
% Accepts:
% 1) paramSetIdx (index in larger sample matrix to start selection of sub-set)
% Saves:
% 1) File: negSumLogL_calib_[date].dat (negative log likelihood for each parameter set in sub-set)

function calib2_sumll4sets(paramSetIdx , tstep_abc , date_abc)

%delete(gcp('nocreate'));

t_curr = tstep_abc;
date = date_abc;

%createParamSet_16Apr20(tstep_abc , date_abc);

%% Cluster information
pc = parcluster('local');    % create a local cluster object
pc.JobStorageLocation = strcat('/gscratch/csde/willmin' , '/' , getenv('SLURM_JOB_ID'))    % explicitly set the JobStorageLocation to the temp directory that was created in the sbatch script
numCPUperNode = str2num(getenv('SLURM_CPUS_ON_NODE'))
parpool(pc , numCPUperNode)    % start the pool with max number workers

%% Load parameters
paramDir = [pwd ,'/Params/'];
paramSetMatrix = load([paramDir,'paramSets_calib_' , date , '_' , num2str(t_curr) , '.dat']);
nPrlSets = 28; %numCPUperNode; %16;
subMatrixInds = [paramSetIdx : (paramSetIdx + nPrlSets - 1)];
pIdx = load([paramDir,'pIdx_calib_' , date , '_0.dat']);

%% If using parameters from a previous calibration or phase, uncomment the following to load resampled subset of parameters from best-fit sets of a previous phase.
%  Note: also need to set paramSet = [ph1sampleSubset(:,subMatrixInds(n)) ; paramSetMatrix(:,subMatrixInds(n))]; or similar, in the correct order of parameters
%  Note: sections to uncomment for Phase 2 in calib1_lhs, calib2_sumll4sets, and abc_smc
pIdx = load([paramDir,'pIdx_calib_' , date , '_0_wPh1Resample.dat']);
ph1sampleSubset = load([paramDir,'resampleSubsetSets_calib_' , date , '_' , num2str(t_curr) , '.dat']);

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
negSumLogLSet = zeros(nPrlSets,1);
parfor n = 1 : nPrlSets
    paramSet = [paramSetMatrix(1:27,subMatrixInds(n)) ; ...
        ph1sampleSubset(1:19,subMatrixInds(n)) ; ...
        paramSetMatrix(28,subMatrixInds(n)) ; ...
        ph1sampleSubset(20:22,subMatrixInds(n)) ; ...
        paramSetMatrix(29:34,subMatrixInds(n))];
    %futureSim(1 , pIdx , paramsSub , paramSet , (paramSetIdx + n - 1) , tstep_abc , date_abc);
    [negSumLogL] = historicalSim(1 , pIdx , paramsSub , paramSet , (paramSetIdx + n - 1) , tstep_abc , date_abc);
    negSumLogLSet(n,1) = negSumLogL;
end

%% Save parameter sets and negSumLogL values
formatOutput = zeros(nPrlSets+(nPrlSets/4) , 1);
paramSetIdxBY4 = [paramSetIdx : 4 : (paramSetIdx + nPrlSets - 1)];
nPrlSetsBY4 = [1 : 4 : (nPrlSets - 1)];
for i = 1 : (nPrlSets/4)
    formatOutput(nPrlSetsBY4(i) + (i-1) : nPrlSetsBY4(i) + 4 + (i-1) , 1) = [paramSetIdxBY4(i) ; negSumLogLSet(((i-1)*4+1 : i*4) , 1)];
end

file = ['negSumLogL_calib_' , date , '_' , num2str(t_curr) , '.dat'];
paramDir = [pwd , '/Params/'];
dlmwrite([paramDir, file] , formatOutput , 'delimiter' , ',' , 'roffset' , 1 , 'coffset' , 0 , '-append' , 'precision' , 9)

%%for j = 1 : nPrlSets
%%    pathModifier = ['toNow_' , date , '_noBaseVax_baseScreen_hpvHIVcalib_' , num2str(t_curr) , '_' , num2str(paramSetIdx + j - 1)];
%%    savdir = [pwd , '/HHCoM_Results/'];
%%    if negSumLogLSet(j,1) < -1000 %-200000.00
%%        delete([savdir , pathModifier , '.mat']);
%%    end
%%end
