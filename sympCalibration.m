sympParams = []; 

for l = 0.0001 : 0.00002 : 0.0002
    for r = 0.014 : 0.002 : 0.02
        for d = 0.6 : 0.1 : 0.8 
            for l_r = 0.02 : 0.02 : 0.1
                for r_d = 0.025 : 0.02 : 0.1

                    sympParams = [sympParams; l r d l_r r_d]; 

                end
            end
        end 
    end 
end 

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

%% Obtain model output for each set of sampled parameters
%
%negSumLogLSet = zeros(nPrlSets,1);
% parfor n = 1 : nPrlSets  
    n=1;

    paramSet = paramSetMatrix(:,subMatrixInds(n));
    %futureSim(1 , pIdx , paramsSub , paramSet , (paramSetIdx + n - 1) , tstep_abc , date_abc);
    

%% Run for loop for the 300 possible symp rate parameter values

parfor sympRun = 1 : length(sympParams)
    sympParams_in = sympParams(sympRun, 1:end);
%     historicalSim(1 , pIdx , paramsSub , paramSet , (paramSetIdx + n - 1) , tstep_abc , date_abc, sympParams_in, sympRun);
    modifiedhistoricalSim(1 , pIdx , paramsSub , paramSet , (paramSetIdx + n - 1) , tstep_abc , date_abc, sympParams_in, sympRun);
end 