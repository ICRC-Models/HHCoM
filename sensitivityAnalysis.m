% Run a global sensitivity analysis
%function sensitivityAnalysis()

%close all; 
clear all; clc
profile clear
profile on

%% Load parameters
paramDir = [pwd ,'\Params\'];
load([paramDir,'settings'])
load([paramDir,'popData'])
load([paramDir,'HIVParams'])
load([paramDir,'general'])
load([paramDir,'mixInfectParams'])
load([paramDir,'vlBeta'])
load([paramDir,'hpvData'])
load([paramDir,'cost_weights'])
load([paramDir,'calibData'])
load([paramDir,'mixInfectIndices'])
load([paramDir,'hivIndices'])
load([paramDir,'hpvIndices'])
load([paramDir,'hpvTreatIndices'])
load([paramDir,'ageRiskInds'])
load([paramDir,'vaxInds'])
load([paramDir,'ager'])
load([paramDir,'vlAdvancer'])
load([paramDir,'fertMat'])
load([paramDir,'hivFertMats'])
load([paramDir,'fertMat2'])
load([paramDir,'hivFertMats2'])
load([paramDir,'vaxer'])
load([paramDir,'circMat'])
load([paramDir,'deathMat'])

%% Latin hypercube sampling of parameter space
nSets = 10;    % number of parameter sets to sample
p = 315;    % number of parameters
sampleNorm = lhsdesign(nSets , p , 'smooth' , 'off');    % latin hypercube sampling

%% Rescale sample values to correct parameter ranges and apply bounds
sampleNorm = sampleNorm';

lb = zeros(p,nSets);
ub = ones(p,nSets);
ub(1:96,:) = 200;
lb(97:99,:) = 0.01;
ub(97:99,:) = 0.1;
lb(100:113,:) = 0.5;
lb(142:155,:) = 0.5;
% ub(114:127,:) = 1.0 - sampleNorm(100:113,:);
% ub(156:169,:) = 1.0 - sampleNorm(142:155,:);
ub(184:279,:) = 200;
lb(280:315,:) = 0.001;
ub(280:291,:) = 0.09;
ub(292:297,:) = 0.5;
ub(298:309,:) = 0.09;
ub(310:315,:) = 0.5;
% KEY
%(1:48):     maleActs, [age,risk], (0.0 to 365) 
%(49:96):    femaleActs, [age,risk], (0.0 to 365)
%(97):       perPartnerHpv, [1x1], (0.0 to 1.0)
%(98):       perPartnerHpv_lr, [1x1], (0.0 to 1.0)
%(99):       perPartnerHpv_nonV, [1x1], (0.0 to 1.0)
%(100:141):  riskDistM, [3:age,risk], (0.0 to 1.0)
%(142:183):  riskDistF, [3:age,risk], (0.0 to 1.0)
%(184:231):  partnersM, [age,risk], (0.0 to 200)
%(232:279):  partnersF, [age,risk], (0.0 to 200)
%(280:297):  muCC, [6x3], (0.0001 to 0.5)
%(298:315):  muCC_det, [6x3], (0.0001 to 0.5)

sample = lb + (sampleNorm .* (ub-lb));
sample(114:127,:) = (1.0 - sample(100:113,:))./2.0 + sampleNorm(114:127,:) .* ((1.0 - sample(100:113,:)) - ((1.0 - sample(100:113,:))./2.0));
sample(156:169,:) = (1.0 - sample(142:155,:))./2.0 + sampleNorm(156:169,:) .* ((1.0 - sample(142:155,:)) - ((1.0 - sample(142:155,:))./2.0));
sample(128:141,:) = 1.0 - sample(100:113,:) - sample(114:127,:);
sample(170:183,:) = 1.0 - sample(142:155,:) - sample(156:169,:);

%% Apply parameter constraints
% indsC1 = any(([sample(132:147,:)<sample(116:131,:)<sample(100:115,:)] < 1) , 1) .* [1:1:nSets];
% sample(:,indsC1) = [];
% if any(size(sample,1)) == 0
%     disp('No samples')
%     error = 1;
% end

%% Obtain model output for each set of sampled parameters
ccIncSet = zeros(nSets,1);
negSumLogLSet = zeros(nSets,1);
parfor n = 1 : nSets
    paramSet = sample(:,n);
    [negSumLogL , ccInc] = calibratorAll(paramSet);
    negSumLogLSet(n,1) = negSumLogL;
    ccIncSet(n,1) = ccInc;
end

figure
subplot(1,2,1);
plot(negSumLogLSet);
title('negSumLogL');
subplot(1,2,2);
plot(ccIncSet);
title('CC Incidence');
