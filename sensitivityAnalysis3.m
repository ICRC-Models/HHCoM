% Generates LHS sample of specified parameters for calibration
% Accepts:
% 1) nSets (number of potential parameter sets to sample)
% Saves:
% 1) File: pIdx_calib_[date].dat (indices for parameters included in
% calibration)
% 2) File: paramSets_calib_[date].dat (sample of parameter sets [number
% parameters x number samples])

function sensitivityAnalysis3(nSets)

%delete(gcp('nocreate')); 
%loadUp(6);

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

pIdx = [1,2,5,6,7,8,9,10,19];    % indices in paramsAll cell array
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
file = 'pIdx_calib_20June19_0.dat';
paramDir = [pwd , '/Params/'];
csvwrite([paramDir, file] , pIdx)

file = 'paramSets_calib_20June19_0.dat';
paramDir = [pwd , '/Params/'];
csvwrite([paramDir, file] , sample)


%%
% KEY
%(1:42):     partnersM, [3:age x risk], (0.001 to 180), *or mults*
%(43:84):    partnersF, [3:age x risk], (0.001 to 180), *or mults*
%(85:126):   riskDistM, [3:age x risk], (0 to 1)
%(127:168):  riskDistF, [3:age x risk], (0 to 1)
%(169):      condUse, [1 x 1], (0.01 to 0.9)
%(170:172):  epsA, [1 x 3], (0.1 to 1)
%(173:175):  epsR, [1 x 3], (0.1 to 1)
%(176:217):  maleActs, [3:age x risk], (1 to 365) --> 25-29 peak acts, *or mults*
%(218:259):  femaleActs, [3:age x risk], (1 to 365) --> 25-29 peak acts, *or mults*
%(260):      perPartnerHpv, [1 x 1], (0.001 to 1.0)
%(261):      perPartnerHpv_lr, [1 x 1], (0.001 to 1.0)
%(262):      perPartnerHpv_nonV, [1 x 1], (0.001 to 1.0)
%(263:266):  hpv_hivMult, [dec CD4 x --hrHPV type--], (0.25x to 4x) 
%(267:282):  rNormal_Inf, [age x 1], (0.25x to 4x) 
%(283:286):  hpv_hivClear, [dec CD4], (0.25x to 4x)
%(287:290):  c3c2Mults, [dec CD4], (0.25x to 4x) 
%(291:294):  c2c1Mults, [dec CD4], (0.25x to 4x)
%(295:297):  kCCDet, [local, reg, dist], (0.001 to 0.5)
%(298:313):  lambdaMultImm, [age x 1], (0.001 to 1)
%(314:315):  maxRateM_vec, [1 x 2], (0.2 to 0.7)
%(316:317):  maxRateF_vec, [1 x 2], (0.2 to 0.7)
%(318):      artHpvMult, [1 x 1], (0.25x to 4x)
%(319:358):  kCD4, [g x vl x CD4], (0.01 to 10) --> ??constraints
%(359:398):  kVL, [g x vl x CD4], (0.01 to 10) --> ??constraints

% lb = ones(p,nSets).*0.001;
% ub = ones(p,nSets);
% ub(1:84,:) = 180;
% lb(85:98,:) = 0.5;
% lb(127:140,:) = 0.5;
% lb(169,:) = 0.01;
% ub(169,:) = 0.9;
% lb(170:175) = 0.1;
% lb(176:259) = 1;
% ub(176:259) = 365;
% lb(263:294) = 0.25;
% ub(263:294) = 4;
% ub(295:297) = 0.5;
% lb(314:315) = 0.2;
% ub(314:315) = 0.7;
% lb(318) = 0.25;
% ub(318) = 4;
% ub(319:398) = 10;

% indsC1 = any(([sample(132:147,:)<sample(116:131,:)<sample(100:115,:)] < 1) , 1) .* [1:1:nSets];
% sample(:,indsC1) = [];
% if any(size(sample,1)) == 0
%     disp('No samples')
%     error = 1;
% end

%% Obtain model output for each set of sampled parameters
% ccIncSet = zeros(nSets,1);
% negSumLogLSet = zeros(nSets,1);
% parfor n = 1 : nSets
%     paramSet = sample(:,n);
%     %[negSumLogL , ccInc] = calibratorAll3(paramSet);
%     [negSumLogL] = calibratorAll3(paramSet);
%     negSumLogLSet(n,1) = negSumLogL;
%     %ccIncSet(n,1) = ccInc;
% end
% 
% figure
% %subplot(1,2,1);
% plot(negSumLogLSet,'o');
% title('negSumLogL');
% % subplot(1,2,2);
% % plot(ccIncSet);
% % title('CC Incidence');

%% Apply parameter constraints
% % riskDistM, riskDistF
% if any(3 == pIdx)
%     idx = find(3 == pIdx);
%     rowL = paramsSub{idx}.length/3;
%     rl = paramsSub{idx}.inds(1:rowL);
%     rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
%     rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
%     sample(rm,:) = (1.0 - sample(rl,:))./2.0 + sampleNorm(rm,:) .* ...
%         ((1.0 - sample(rl,:)) - ((1.0 - sample(rl,:))./2.0));
%     sample(rh,:) = 1.0 - sample(rl,:) - sample(rm,:);
% end
% if any(4 == pIdx)
%     idx = find(4 == pIdx);
%     rowL = paramsSub{idx}.length/3;
%     rl = paramsSub{idx}.inds(1:rowL);
%     rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
%     rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
%     sample(rm,:) = (1.0 - sample(rl,:))./2.0 + sampleNorm(rm,:) .* ...
%         ((1.0 - sample(rl,:)) - ((1.0 - sample(rl,:))./2.0));
%     sample(rh,:) = 1.0 - sample(rl,:) - sample(rm,:);
% end

%%
% file = 'negSumLogL_calib_25Feb19.dat';
% paramDir = [pwd , '\Params\'];
%csvwrite([paramDir, file] , negSumLogLSet)
