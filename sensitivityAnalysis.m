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
nSets = 32;
p = 241;
sampleNorm = lhsdesign(nSets , p , 'smooth' , 'off');

%% Rescale sample values to correct parameter ranges
lb = zeros(241,1) * 0.0001;
ub = ones(241,1);
ub(9:104) = 200.0;
lb(112:119) = 0.5;
ub(112:119) = 4.0;
lb(128:135) = 1.0;
ub(128:131) = c3c2Mults .* 4;
ub(132:135) = c2c1Mults .* 4;
%ub(157:236) = 10.0;
lb(241) = 1.0;
ub(241) = 4.0;

sampleNorm = sampleNorm';
sample = lb + (sampleNorm .* (ub-lb));

%% Obtain model output for each set of sampled parameters
ccIncSet = zeros(nSets,1);
negSumLogLSet = zeros(nSets,1);
parfor n = 1 : nSets
    paramSet = sample(:,n);
    [negSumLogL , ccInc] = calibratorAll(paramSet);
    negSumLogL(n,1) = negSumLogL;
    ccIncSet(n,1) = ccInc;
end

figure
subplot(1,2,1);
plot(negSumLogLSet);
subplot(1,2,2);
plot(ccIncSet);
