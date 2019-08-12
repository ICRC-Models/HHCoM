% Identify parameter ranges of top fifty best-fitting sets and compare to bounds
function [] = idParamRanges(tstep_abc , date_abc)
t_curr = tstep_abc;
date = date_abc;

%% Load all particles
paramDir = [pwd , '/Params/'];
masterSetMatrix = load([paramDir , 'masterSets_calib_' , date , '_' , num2str(t_curr) , '.dat']); % load most recent parameter sample
pIdx = load([paramDir , 'pIdx_calib_' , date , '_0' , '.dat']); % load parameter indices
orderedLL = load([paramDir , 'orderedLL_calib_' , date , '_' , num2str(t_curr) , '.dat']); % load most recent ordered log-likelihoods

%% Get indices and parameter values of top fifty sets
top50Inds = orderedLL(1:50,1);
top50Params = masterSetMatrix(:,top50Inds);

%% Find min and max of each parameter and compare to bounds
[paramsAll] = genParamStruct();
paramsSub = cell(length(pIdx),1);
lb = [];
ub = [];
for s = 1 : length(pIdx)
    paramsSub{s} = paramsAll{pIdx(s)};
    lb = [lb; paramsSub{s}.lb];
    ub = [ub; paramsSub{s}.ub];
end

minVals = zeros(size(top50Params,1),1);
maxVals = minVals;
lbVals = minVals;
ubVals = minVals;
for i = 1 : size(top50Params,1)
    minVals(i,1) = min(top50Params(i,:));
    maxVals(i,1) = max(top50Params(i,:));
    lbVals(i,1) = lb(i);
    ubVals(i,1) = ub(i);
end

%% Print ranges
fileRanges = ['rangesVSbounds_calib_' , date , '_' , num2str(t_curr) , '.dat']; % save ranges compared to bounds
csvwrite([paramDir, fileRanges , [minVals , maxVals , lbVals, ubVals]);


