% Parameter estimation and sensitivity analysis
function sensitivityAnalysis3_prt2(paramSetIdx)

%% Load parameters
paramDir = [pwd ,'/Params/'];
load([paramDir,'settings']);
load([paramDir,'general']);
paramSetMatrix = load([paramDir,'paramSets_calib_02May19.dat']);
nPrlSets = 16;
subMatrixInds = [paramSetIdx : (paramSetIdx + nPrlSets)];

%% Cluster information
numCores = feature('numcores');
parpool(numCores);

%% Obtain model output for each set of sampled parameters
%ccIncSet = zeros(nSets,1);
negSumLogLSet = zeros(nPrlSets,1);
parfor n = 1 : nPrlSets
    paramSet = paramSetMatrix(:,subMatrixInds(n));
    [negSumLogL] = calibratorAll3(paramSet);
    negSumLogLSet(n,1) = negSumLogL;
end

%% Save parameter sets and negSumLogL values
file = 'negSumLogL_calib_02May19.dat';
paramDir = [pwd , '/Params/'];
dlmwrite([paramDir, file] , [paramSetIdx; negSumLogLSet] , 'delimiter' , ',' , 'roffset' , 1 , 'coffset' , 0 , '-append')
