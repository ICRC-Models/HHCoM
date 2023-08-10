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

%% Select the kSymp parameters

file = [pwd , '/Params/kSymp_ParamSets_recalib_round2.xlsx'];

paramSet_kSymp = xlsread(file, 'Param sets', 'A2:C101'); 

%% Obtain model output for each set of sampled parameters
parfor sympRun = 1 : length(paramSet_kSymp)
    sympParams_in = paramSet_kSymp(sympRun, 1:end);

    m = 1; % originally it was paramSetIdx:nPrlSets:(numBestFits+paramSetIdx-1), but this just is 1
    subMatrixInds = [m : (m + nPrlSets - 1)];
    
    % select a random parameter set from 1 to 25
    n = randi([1, 25]);
    paramSet = top50Params(:,subMatrixInds(n));

    modifiedhistoricalSim(1 , pIdx , paramsSub , paramSet , (paramSetIdx + n - 1) , tstep_abc , date_abc, sympParams_in, sympRun);

end 
