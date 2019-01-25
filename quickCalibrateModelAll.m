function quickCalibrateModel()

close all; clear all; clc

%% Load parameters
paramDir = [pwd ,'\Params\'];
load([paramDir,'settings'])
load([paramDir,'popData'])
load([paramDir,'HIVParams'])
load([paramDir,'general'])
load([paramDir,'mixInfectParams'])
load([paramDir,'vlBeta'])
load([paramDir,'hpvData'])
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

import java.util.LinkedList

%% Set parameter upper and lower bounds 

lb = zeros(345,1);
ub = ones(345,1);

ub(9:104) = 365.0;
lb(112:127) = kCin1_Inf ./ 10;
lb(128:143) = kCin2_Cin1 ./ 10;
lb(144:159) = kCin3_Cin2 ./ 10;
lb(160:175) = kCC_Cin3 ./ 10;
lb(176:191) = rNormal_Inf ./ 10;
lb(192:207) = kInf_Cin1 ./ 10;
lb(208:223) = kCin1_Cin2 ./ 10;
lb(224:239) = kCin2_Cin3 ./ 10;
ub(112:127) = kCin1_Inf .* 10;
ub(128:143) = kCin2_Cin1 .* 10;
ub(144:159) = kCin3_Cin2 .* 10;
ub(160:175) = kCC_Cin3 .* 10;
ub(176:191) = rNormal_Inf .* 10;
ub(192:207) = kInf_Cin1 .* 10;
ub(208:223) = kCin1_Cin2 .* 10;
ub(224:239) = kCin2_Cin3 .* 10;
lb(248:255) = 1.0;
ub(248:251) = c3c2Mults .* 10;
ub(252:255) = c2c1Mults .* 10;
ub(277:340) = 10.0;
lb(345) = 1.0;
ub(345) = 10.0;

%(1:3):     epsA, [3x1], (0.0 to 1.0), reset in mixInfect
%(4:6):     epsR, [3x1], (0.0 to 1.0), reset in mixInfect
%(7):       prepOut, [1x1], (0.0 to 1.0), reset in hiv2a to 0
%(8):       artOut, [1x1], (0.0 to 1.0), reset in hiv2a to 0
%(9:56):    maleActs, [age,risk], (0.0 to 365) 
%(57:104):  femaleActs, [age,risk], (0.0 to 365)
%(105):     perPartnerHpv, [1x1], (0.0 to 1.0)
%(106):     perPartnerHpv_lr, val, (0.0 to 1.0)
%(107):     perPartnerHpv_nonV, val, (0.0 to 1.0)
%(108:111): hpv_hivMult, [CD4x1], (0.0 to 1.0) 
%(112:127): kCin1_Inf, [agex1], init, (/10 to x10) or (0.0 to 1.0)
%(128:143): kCin2_Cin1, [agex1], init, (/10 to x10)
%(144:159): kCin3_Cin2, [agex1], init, (/10 to x10)
%(160:175): kCC_Cin3, [agex1], init, (/10 to x10)
%(176:191): rNormal_Inf, [agex1], init, (/10 to x10)
%(192:207): kInf_Cin1, [agex1], init, (/10 to x10)
%(208:223): kCin1_Cin2, [agex1], init, (/10 to x10)
%(224:239): kCin2_Cin3, [agex1], init, (/10 to x10)
%(240:243): hpv_hivClear, [CD4x1], (0.0 to 1.0)
%(244:247): rImmuneHiv, [CD4x1], (0.0 to 1.0)
%(248:251): c3c2Mults, [CD4x1], init, (all ones to x10)
%(252:255): c2c1Mults, [CD4x1], init, (all ones to x10)
%(256:271): lambdaMultImm, [agex1], (0.0 to 1.0)
%(272):     kRL, [1x1], (0.0 to 1.0)
%(273):     kDR, [1x1], (0.0 to 1.0)
%(274:276): kCCDet, [3x1], (0.0 to 1.0)
%(277:308): kCD4, [genderxvlxcd4], (0.0 to 10.0)
%(309:340): kVL, [genderxcd4xvl], (0.0 to 10.0)
%(341:342): maxRateM_vec,[2x1], (0.0 to 1.0), reset in mainCalibrated
%(343:344): maxRateF_vec, [2x1], (0.0 to 1.0), reset in mainCalibrated
%(345):     artHpvMult, [1x1], (1.0 to x10)

%% Set parameter constraints
% max HPV acquisition multiplier for HIV+ on ART <= max acquisition multiplier for HIV+
% restrict CD4 progressions to be highest for the lowest CD4 category because 
%   MLE might suggest progressions that were not in the same rank order of CD4 categories

%% Find optimal parameter set
% loop over multiple start points, then compare local optima to find global
% optimum

numStartPts = 10;
x = zeros(length(lb),numStartPts);

for i = 1 : numStartPts
    initParams = lb + rand(size(lb)).*(ub - lb);

    options = psoptimset('UseParallel' , true , 'Cache' , 'on' ,...
        'CacheTol' , 0.1 , 'CompletePoll' , 'on' , 'TolMesh' , 0.1, ...
        'Display','iter','PlotFcn',@psplotbestf);
    x(:,i) = patternsearch(@calibratorAll, initParams , [] , [] , [] , [] , lb , ub , [] , options);
end
%% Save calibrated parameters
file = 'HPV_calib_23Jan19.dat';
paramDir = [pwd , '\Params\'];
csvwrite([paramDir, file] , x)
