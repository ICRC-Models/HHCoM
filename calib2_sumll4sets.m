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

%createParamSet(tstep_abc , date_abc);

%% Cluster information
pc = parcluster('local');    % create a local cluster object
pc.JobStorageLocation = strcat('/gscratch/csde/carajb' , '/' , getenv('SLURM_JOB_ID'))    % explicitly set the JobStorageLocation to the temp directory that was created in the sbatch script
numCPUperNode = str2num(getenv('SLURM_CPUS_ON_NODE'))
parpool(pc , numCPUperNode)    % start the pool with max number workers

%% Load parameters
paramDir = [pwd ,'/Params/'];
paramSetMatrix = load([paramDir,'paramSets_calib_' , date , '_' , num2str(t_curr) , '_7867mod45' , '.dat']);
nPrlSets = 1; %28; %numCPUperNode; %16;
subMatrixInds = [paramSetIdx : (paramSetIdx + nPrlSets - 1)];
pIdx = load([paramDir,'pIdx_calib_' , date , '_0.dat']);

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
for n = 1 : nPrlSets
    paramSet = paramSetMatrix(:,subMatrixInds(n));
    %futureSim(1 , pIdx , paramsSub , paramSet , (paramSetIdx + n - 1) , tstep_abc , date_abc);
    [negSumLogL] = historicalSim(1 , pIdx , paramsSub , paramSet , (paramSetIdx + n - 1) , tstep_abc , date_abc);
    negSumLogLSet(n,1) = negSumLogL
end

%% Save parameter sets and negSumLogL values
% formatOutput = zeros(nPrlSets+(nPrlSets/4) , 1);
% paramSetIdxBY4 = [paramSetIdx : 4 : (paramSetIdx + nPrlSets - 1)];
% nPrlSetsBY4 = [1 : 4 : (nPrlSets - 1)];
% for i = 1 : (nPrlSets/4)
%     formatOutput(nPrlSetsBY4(i) + (i-1) : nPrlSetsBY4(i) + 4 + (i-1) , 1) = [paramSetIdxBY4(i) ; negSumLogLSet(((i-1)*4+1 : i*4) , 1)];
% end
% 
% file = ['negSumLogL_calib_' , date , '_' , num2str(t_curr) , '.dat'];
% paramDir = [pwd , '/Params/'];
% dlmwrite([paramDir, file] , formatOutput , 'delimiter' , ',' , 'roffset' , 1 , 'coffset' , 0 , '-append' , 'precision' , 9)
