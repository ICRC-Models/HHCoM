% Generates LHS sample of specified parameters for calibration
% Accepts:
% 1) nSets (number of potential parameter sets to sample)
% Saves:
% 1) File: pIdx_calib_[date].dat (indices for parameters included in
% calibration)
% 2) File: paramSets_calib_[date].dat (sample of parameter sets [number
% parameters x number samples])

function calib1_lhs(nSets , tstep_abc , date_abc)

%delete(gcp('nocreate')); 

t_curr = tstep_abc;
date = date_abc;

%% Cluster information
pc = parcluster('local');    % create a local cluster object
pc.JobStorageLocation = strcat('/gscratch/csde/willmin' , '/' , getenv('SLURM_JOB_ID'))    % explicitly set the JobStorageLocation to the temp directory that was created in the sbatch script
parpool(pc , str2num(getenv('SLURM_CPUS_ON_NODE')))    % start the pool with max number workers

%% Load structure of all potentially calibrated parameters
[paramsAll] = genParamStruct();

%% Latin hypercube sampling of parameter space
pIdx = [1,2,5,6,7,8,9,10,35,38];    % indices in paramsAll cell array

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

sampleNorm = lhsdesign(nSets , p , 'smooth' , 'off');    % latin hypercube sampling

%% Rescale sample values to correct parameter ranges and apply bounds
sampleNorm = sampleNorm';

sample = lb + (sampleNorm .* (ub-lb));

%% Save parameter sets and parameter index values
file = ['pIdx_calib_' , date , '_' , num2str(t_curr) , '.dat'];
paramDir = [pwd , '/Params/'];
csvwrite([paramDir, file] , pIdx)

file = ['paramSets_calib_' , date , '_' , num2str(t_curr) , '.dat'];
paramDir = [pwd , '/Params/'];
csvwrite([paramDir, file] , sample)

%% If on Phase 2 of calibration, uncomment the following to resample a subset of parameters from best-fit sets of a previous phase.
%  Note: sections to uncomment for Phase 2 in calib1_lhs, calib2_sumll4sets, and abc_smc
% pIdx_wPh1Resample = [1,2,5,6,9,35, pIdx(1,:)];
% 
% ph1_top50Sets = load([paramDir,'alphaParamSets_calib_22Apr20_20_top50Sets.dat']);
% ph1sample = datasample(ph1_top50Sets, nSets , 2); % resample
% ph1sampleSubset = [ph1sample(1:21,:); ph1sample(26,:)]; % keep subset of resampled parameter set
% 
% file = ['pIdx_calib_' , date , '_' , num2str(t_curr) , '_wPh1Resample' , '.dat'];
% paramDir = [pwd , '/Params/'];
% csvwrite([paramDir, file] , pIdx_wPh1Resample)
% 
% file = ['resampleSets_calib_' , date , '_' , num2str(t_curr) , '.dat'];
% paramDir = [pwd , '/Params/'];
% csvwrite([paramDir, file] , ph1sample)
% 
% file = ['resampleSubsetSets_calib_' , date , '_' , num2str(t_curr) , '.dat'];
% paramDir = [pwd , '/Params/'];
% csvwrite([paramDir, file] , ph1sampleSubset)
