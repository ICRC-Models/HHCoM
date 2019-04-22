% Parameter estimation and sensitivity analysis
%function sensitivityAnalysis()

close all; clear all; clc

%loadUp(6);

%% Load parameters
paramDir = [pwd ,'\Params\'];
load([paramDir,'settings'])
load([paramDir,'general'])

%% Cluster information
% Create a local cluster object
myCluster = parcluster('local'); 
% Set the JobStorageLocation to the temporary directory that was created in your slurm script
myCluster.JobStorageLocation = strcat('/gscratch/csde/carajb/', getenv('SLURM_JOB_ID')) ;
numCores = feature('numcores')
parpool(numCores);


%% Latin hypercube sampling of parameter space
nSets = 48; %100;    % number of parameter sets to sample
p = 84; %398;    % number of parameters
sampleNorm = lhsdesign(nSets , p , 'smooth' , 'off');    % latin hypercube sampling

%% Rescale sample values to correct parameter ranges and apply bounds
sampleNorm = sampleNorm';

lb = ones(p,nSets).*0.001;
ub = ones(p,nSets);
ub(1:84,:) = 180;
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

% KEY
%(1:42):     partnersM, [3:age x risk], (0.001 to 180)
%(43:84):    partnersF, [3:age x risk], (0.001 to 180)
%(85:126):   riskDistM, [3:age x risk], (0 to 1)
%(127:168):  riskDistF, [3:age x risk], (0 to 1)
%(169):      condUse, [1 x 1], (0.01 to 0.9)
%(170:172):  epsA, [1 x 3], (0.1 to 1)
%(173:175):  epsR, [1 x 3], (0.1 to 1)
%(176:217):  maleActs, [3:age x risk], (1 to 365) --> 25-29 peak acts
%(218:259):  femaleActs, [3:age x risk], (1 to 365) --> 25-29 peak acts
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

sample = lb + (sampleNorm .* (ub-lb));

%% Apply parameter constraints
% partnersM, partnersF
sample(15:28,:) = (sample(29:42,:)-lb(15:28,:))./2.0 + sampleNorm(15:28,:) .* (sample(29:42,:) - ((sample(29:42,:)-lb(15:28,:))./2.0));
sample(57:70,:) = (sample(71:84,:)-lb(57:70,:))./2.0 + sampleNorm(57:70,:) .* (sample(71:84,:) - ((sample(71:84,:)-lb(57:70,:))./2.0));
sample(1:14,:) = lb(1:14,:) + sampleNorm(1:14,:) .* (sample(15:28,:) - lb(1:14,:));
sample(43:56,:) = lb(43:56,:) + sampleNorm(43:56,:) .* (sample(57:70,:) - lb(43:56,:));
% % riskDistM, riskDistF
% sample(99:112,:) = (1.0 - sample(85:98,:))./2.0 + sampleNorm(99:112,:) .* ((1.0 - sample(85:98,:)) - ((1.0 - sample(85:98,:))./2.0));
% sample(141:154,:) = (1.0 - sample(127:140,:))./2.0 + sampleNorm(141:154,:) .* ((1.0 - sample(127:140,:)) - ((1.0 - sample(127:140,:))./2.0));
% sample(113:126,:) = 1.0 - sample(85:98,:) - sample(99:112,:);
% sample(155:168,:) = 1.0 - sample(127:140,:) - sample(141:154,:);
% % maleActs, femaleActs
% sample(190:203,:) = (sample(176:189,:)-lb(190:203,:))./2.0 + sampleNorm(190:203,:) .* (sample(176:189,:) - ((sample(176:189,:)-lb(190:203,:))./2.0));
% sample(204:217,:) = lb(204:217,:) + sampleNorm(204:217,:) .* (sample(190:203,:) - lb(204:217,:));
% sample(232:245,:) = (sample(218:231,:)-lb(232:245,:))./2.0 + sampleNorm(232:245,:) .* (sample(218:231,:) - ((sample(218:231,:)-lb(232:245,:))./2.0));
% sample(246:259,:) = lb(246:259,:) + sampleNorm(246:259,:) .* (sample(232:245,:) - lb(246:259,:));

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

%% Save parameter sets and negSumLogL values
file = 'paramSets_calib_22Apr19.dat';
paramDir = [pwd , '\Params\'];
csvwrite([paramDir, file] , sample)
% file = 'negSumLogL_calib_25Feb19.dat';
% paramDir = [pwd , '\Params\'];
% csvwrite([paramDir, file] , negSumLogLSet)
