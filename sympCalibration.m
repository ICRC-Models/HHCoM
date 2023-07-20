paramSetIdx=1;
tstep_abc=6;
date_abc='22Apr20Ph2V11';
username='clh89';

t_curr = tstep_abc;
date = date_abc;

%% Cluster information -- Erisone
pc = parcluster('local'); 
pc.JobStorageLocation = getenv('TMPDIR'); % how to pull job id?
numCPUperNode = 28; % how to pull CPUs on node? set to 8 as an initial test.
parpool(pc, numCPUperNode)

%% Load parameters
paramDir = [pwd ,'/Params/'];
paramSetMatrix = load([paramDir,'stochasticParamsets.dat']);
nPrlSets = 25; %numCPUperNode; %16;
subMatrixInds = [paramSetIdx : (paramSetIdx + nPrlSets - 1)];
pIdx = [10, 35, 38];

[paramsAll] = genParamStruct();
paramsSub = cell(length(pIdx),1);
startIdx = 1;
for s = 1 : length(pIdx)
    paramsSub{s}.length = paramsAll{pIdx(s)}.length;
    paramsSub{s}.inds = (startIdx : (startIdx + paramsSub{s}.length - 1));
    startIdx = startIdx + paramsSub{s}.length;
end

%% Select the kSymp parameters

file = [pwd , '/Params/kSymp_ParamSets_recalib_round2.xlsx'];

paramSet_kSymp = xlsread(file, 'Param sets', 'A2:C101'); 
    

%% Run parfor loop for all the parameter sets 

parfor sympRun = 1 : length(paramSet_kSymp)
    sympParams_in = paramSet_kSymp(sympRun, 1:end);

    % select a random parameter from paramsSub
    n = randi([1, 25]); 
    paramSet = paramSetMatrix(:,subMatrixInds(n));

    modifiedhistoricalSim(1 , pIdx , paramsSub , paramSet , (paramSetIdx + n - 1) , tstep_abc , date_abc, sympParams_in, sympRun);
end 