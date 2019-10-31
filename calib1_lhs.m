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
%loadUp(6);

t_curr = tstep_abc;
date = date_abc;

%% Load parameters
paramDir = [pwd ,'/Params/'];
load([paramDir,'settings'])
load([paramDir,'general'])

%% Cluster information
pc = parcluster('local');    % create a local cluster object
pc.JobStorageLocation = strcat('/gscratch/csde/carajb' , '/' , getenv('SLURM_JOB_ID'))    % explicitly set the JobStorageLocation to the temp directory that was created in the sbatch script
parpool(pc , str2num(getenv('SLURM_CPUS_ON_NODE')))    % start the pool with max number workers

%% Load structure of all potentially calibrated parameters
[paramsAll] = genParamStruct();

%% Latin hypercube sampling of parameter space
%nSets = 48; %100;    % number of parameter sets to sample
%p = 84; 398;    % number of parameters

pIdx = [1,2,5,6,7,8,9,10,11,18,24];    % indices in paramsAll cell array
prtnrActMults = 1;

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

%% Apply parameter constraints

% partnersM, partnersF (if calibrating actual values vs. multipliers)
if (any(1 == pIdx) && ~prtnrActMults)
    idx = find(1 == pIdx);
    rowL = paramsSub{idx}.length/3;
    rl = paramsSub{idx}.inds(1:rowL);
    rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    sample(rm,:) = (sample(rh,:)-lb(rm,:))./2.0 + sampleNorm(rm,:) .* ...
        (sample(rh,:) - ((sample(rh,:)-lb(rm,:))./2.0)); % mr partners < hr partners, > lr partners
    sample(rl,:) = lb(rl,:) + sampleNorm(rl,:) .* (sample(rm,:) - lb(rl,:)); % lr partners < mr partners
end
if (any(2 == pIdx) && ~prtnrActMults) 
    idx = find(2 == pIdx);
    rowL = paramsSub{idx}.length/3;
    rl = paramsSub{idx}.inds(1:rowL);
    rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    sample(rm,:) = (sample(rh,:)-lb(rm,:))./2.0 + sampleNorm(rm,:) .* ...
        (sample(rh,:) - ((sample(rh,:)-lb(rm,:))./2.0));
    sample(rl,:) = lb(rl,:) + sampleNorm(rl,:) .* (sample(rm,:) - lb(rl,:));
end

% maleActs, femaleActs (if calibrating actual values vs. multipliers)
if (any(8 == pIdx) && ~prtnrActMults)
    idx = find(8 == pIdx);
    rowL = paramsSub{idx}.length/3;
    rl = paramsSub{idx}.inds(1:rowL);
    rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    sample(rm,:) = (sample(rl,:)-lb(rm,:))./2.0 + sampleNorm(rm,:) .* ...
        (sample(rl,:) - ((sample(rl,:)-lb(rm,:))./2.0));
    sample(rh,:) = lb(rh,:) + sampleNorm(rh,:) .* (sample(rm,:) - lb(rh,:));
end
if (any(9 == pIdx) && ~prtnrActMults)
    idx = find(9 == pIdx);
    rowL = paramsSub{idx}.length/3;
    rl = paramsSub{idx}.inds(1:rowL);
    rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    sample(rm,:) = (sample(rl,:)-lb(rm,:))./2.0 + sampleNorm(rm,:) .* ...
        (sample(rl,:) - ((sample(rl,:)-lb(rm,:))./2.0));
    sample(rh,:) = lb(rh,:) + sampleNorm(rh,:) .* (sample(rm,:) - lb(rh,:));
end

%% Save parameter sets and negSumLogL values
file = ['pIdx_calib_' , date , '_' , num2str(t_curr) , '.dat'];
paramDir = [pwd , '/Params/'];
csvwrite([paramDir, file] , pIdx)

file = ['paramSets_calib_' , date , '_' , num2str(t_curr) , '.dat'];
paramDir = [pwd , '/Params/'];
csvwrite([paramDir, file] , sample)
