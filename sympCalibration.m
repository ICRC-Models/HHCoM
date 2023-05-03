sympParams = []; 

for l = 0.0001 : 0.0004 : 0.001
    for r = 0.001 : 0.004 : 0.1
        for d = 0.7 : 0.1 : 1.0 

            sympParams = [sympParams; l r d]; 

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

%%
nPrlSets = 25;

%% Load all particles
paramDir = [pwd , '/Params/'];
% masterSetMatrix = load(string(strjoin([paramDir , 'masterSets_calib_' , date , '_' , num2str(t_curr) , '.dat'],''))); % load most recent parameter sample
masterSetMatrix = load([paramDir , 'masterSets_calib_' , date , '_' , num2str(t_curr) , '.dat']); 
% pIdx = load(string(strjoin([paramDir , 'pIdx_calib_' , date , '_0' , '.dat'],''))); % load parameter indices
pIdx = load([paramDir , 'pIdx_calib_' , date , '_0' , '.dat']);
% orderedLL = load(string(strjoin([paramDir , 'orderedLL_calib_' , date , '_' , num2str(t_curr) , '.dat'],''))); % load most recent ordered log-likelihoods
orderedLL = load([paramDir , 'orderedLL_calib_' , date , '_' , num2str(t_curr) , '.dat']);

%% Get indices and parameter values of numBestFits*2  sets
numBestFits = 25;
numSets50 = size(orderedLL , 1)*((numBestFits*2)/5600);
specIndsList = [1,2,3,6,8,9,11,12,13,15,20,21,22,26,27,32,34,35,38,39,40,41,42,45,47]; 
top50Inds = orderedLL(specIndsList,1); %orderedLL(1:numSets50,1);
top50Params = masterSetMatrix(:,top50Inds);

%% If on Phase 2 of calibration, uncomment the following to plot resampled Ph1 parameters + Ph2 parameters
% pIdx = load(string(strjoin([paramDir,'pIdx_calib_' , date , '_0_wPh1Resample.dat'],'')));
pIdx = load([paramDir,'pIdx_calib_' , date , '_0_wPh1Resample.dat']);
% masterResampleSubsetMatrix = load(string(strjoin([paramDir , 'masterResampleSubsetMatrix_calib_' , date , '_' , num2str(t_curr) , '.dat'],''))); % load most recent Ph1 resampled parameters
masterResampleSubsetMatrix = load([paramDir , 'masterResampleSubsetMatrix_calib_' , date , '_' , num2str(t_curr) , '.dat']);
masterCombinedPhaseMatrix = [masterResampleSubsetMatrix ; masterSetMatrix];
top50Params = masterCombinedPhaseMatrix(:,top50Inds);

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
m=1;
n=1; 
subMatrixInds = [m : (m + nPrlSets - 1)];
paramSet = top50Params(:,subMatrixInds(n));

%% Run for loop for the 300 possible symp rate parameter values

parfor sympRun = 1 : length(sympParams)
    sympParams_in = sympParams(sympRun, :);
    modifiedhistoricalSim(1 , pIdx , paramsSub , paramSet ,specIndsList(m + n - 1) , tstep_abc , date_abc, sympParams_in, sympRun);
end 