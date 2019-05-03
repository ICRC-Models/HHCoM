% Parameter estimation and sensitivity analysis
function sensitivityAnalysis3_prt2(paramSetIdx)

delete(gcp('nocreate')); 

%% Load parameters
paramDir = [pwd ,'/Params/'];
load([paramDir,'settings']);
load([paramDir,'general']);
paramSetMatrix = load([paramDir,'paramSets_calib_02May19.dat']);
nPrlSets = 16;
subMatrixInds = [paramSetIdx : (paramSetIdx + nPrlSets)];
pIdx = load([paramDir,'pIdx_calib_02May19.dat']);

[paramsAll] = genParamStruct();
paramsSub = cell(length(pIdx),1);
startIdx = 1;
for s = 1 : length(pIdx)
    paramsSub{s}.length = paramsAll{s}.length;
    paramsSub{s}.inds = [(startIdx : (startIdx + paramsAll{s}.length - 1)) , 1];
    startIdx = startIdx + paramsAll{s}.length;
end

%% Cluster information
numCores = feature('numcores');
parpool(numCores);

%% Obtain model output for each set of sampled parameters
%ccIncSet = zeros(nSets,1);
negSumLogLSet = zeros(nPrlSets,1);
for n = 1 : nPrlSets
    paramSet = paramSetMatrix(:,subMatrixInds(n));
    [negSumLogL] = calibratorAll3(pIdx , paramsSub , paramSet);
    negSumLogLSet(n,1) = negSumLogL;
end

%% Save parameter sets and negSumLogL values
file = 'negSumLogL_calib_02May19.dat';
paramDir = [pwd , '/Params/'];
csvwrite([paramDir, file] , negSumLogLSet,(paramSetIdx-1),0)
