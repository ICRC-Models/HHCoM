% Generates LHS sample of specified parameters for calibration
% Accepts:
% 1) nSets (number of potential parameter sets to sample)
% Saves:
% 1) File: pIdx_calib_[date].dat (indices for parameters included in
% calibration)
% 2) File: paramSets_calib_[date].dat (sample of parameter sets [number
% parameters x number samples])

function lhsPatternSrch(tstep_abc , date_abc)

%delete(gcp('nocreate')); 

t_curr = tstep_abc;
date = date_abc;

%% Cluster information
pc = parcluster('local');    % create a local cluster object
pc.JobStorageLocation = strcat('/gscratch/csde/carajb' , '/' , getenv('SLURM_JOB_ID'))    % explicitly set the JobStorageLocation to the temp directory that was created in the sbatch script
parpool(pc , str2num(getenv('SLURM_CPUS_ON_NODE')))    % start the pool with max number workers

%% Load structure of all potentially calibrated parameters
[paramsAll] = genParamStruct();

%% Latin hypercube sampling of parameter space
nSets = 1;
pIdx = [1,2,5,6,7,8,9,10,18];    % indices in paramsAll cell array

paramsSub = cell(length(pIdx),1);
p = 0;
startIdx = 1;
lb = [];
ub = [];
for s = 1 : length(pIdx)
    paramsSub{s} = paramsAll{pIdx(s)};
    p = p + paramsSub{s}.length;
    paramsSub{s}.inds = (startIdx : (startIdx + paramsSub{s}.length - 1));
    startIdx = startIdx + paramsSub{s}.length;

    lb = [lb; paramsSub{s}.lb];
    ub = [ub; paramsSub{s}.ub];
end
lb(8) = lb(8).*10; % re-scale perPartnerHpv to be more similar in scale to other params
ub(8) = ub(8).*10;

sampleNorm = lhsdesign(nSets , p , 'smooth' , 'off');    % latin hypercube sampling

%% Rescale sample values to correct parameter ranges and apply bounds
sampleNorm = sampleNorm';

sample = lb + (sampleNorm .* (ub-lb));

%% Save parameter sets
file = ['pIdx_patternSrch_' , date , '_' , num2str(t_curr) , '.dat'];
paramDir = [pwd , '/Params/'];
csvwrite([paramDir, file] , pIdx)

file = ['paramSets_patternSrch_' , date , '_' , num2str(t_curr) , '.dat'];
paramDir = [pwd , '/Params/'];
csvwrite([paramDir, file] , sample)

