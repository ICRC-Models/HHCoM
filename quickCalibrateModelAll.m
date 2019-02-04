function quickCalibrateModel()

close all; clear all; clc
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

%% Set parameter upper and lower bounds 
lb = zeros(241,1);
ub = ones(241,1);
ub(9:104) = 365.0;
lb(112:119) = 0.001;
ub(112:119) = 20.0;
lb(128:135) = 1.0;
ub(128:131) = c3c2Mults .* 10;
ub(132:135) = c2c1Mults .* 10;
%ub(157:236) = 10.0;
lb(241) = 1.0;
ub(241) = 10.0;

% KEY
%(1:3):     epsA, [3x1], (0.0 to 1.0), XXXreset in mixInfectXXX
%(4:6):     epsR, [3x1], (0.0 to 1.0), XXXreset in mixInfectXXX
%(7):       prepOut, [1x1], (0.0 to 1.0), reset in hiv2a to 0
%(8):       artOut, [1x1], (0.0 to 1.0), reset in hiv2a to 0
%(9:56):    maleActs, [age,risk], (0.0 to 365) 
%(57:104):  femaleActs, [age,risk], (0.0 to 365)
%(105):     perPartnerHpv, [1x1], (0.0 to 1.0)
%(106):     perPartnerHpv_lr, [1x1], (0.0 to 1.0)
%(107):     perPartnerHpv_nonV, [1x1], (0.0 to 1.0)
%(108:111): hpv_hivMult, [CD4x1], (0.0 to 1.0) 
%(112): kCin1_Inf, [1x1], (0.001 to 20)
%(113): kCin2_Cin1, [1x1], (0.001 to 20)
%(114): kCin3_Cin2, [1x1], (0.001 to 20)
%(115): kCC_Cin3, [1x1], (0.001 to 20)
%(116): rNormal_Inf, [1x1], (0.001 to 20)
%(117): kInf_Cin1, [1x1], (0.001 to 20)
%(118): kCin1_Cin2, [1x1], (0.001 to 20)
%(119): kCin2_Cin3, [1x1], (0.001 to 20)
%(120:123): hpv_hivClear, [CD4x1], (0.0 to 1.0)
%(124:127): rImmuneHiv, [CD4x1], (0.0 to 1.0)
%(128:131): c3c2Mults, [CD4x1], init, (all ones to x10)
%(132:135): c2c1Mults, [CD4x1], init, (all ones to x10)
%(136:151): lambdaMultImm, [agex1], (0.0 to 1.0)
%(152):     kRL, [1x1], (0.0 to 1.0)
%(153):     kDR, [1x1], (0.0 to 1.0)
%(154:156): kCCDet, [3x1], (0.0 to 1.0)
%(157:196): kCD4, [genderxvlxcd4], (0.0 to 10.0)
%(197:236): kVl, [genderxcd4xvl], (0.0 to 10.0)
%(237:238): maxRateM_vec,[2x1], (0.0 to 1.0), reset in mainCalibrated
%(239:240): maxRateF_vec, [2x1], (0.0 to 1.0), reset in mainCalibrated
%(241):     artHpvMult, [1x1], (1.0 to x10)

%% Set parameter constraints
% max HPV acquisition multiplier for HIV+ on ART <= max acquisition multiplier for HIV+
% restrict CD4 progressions to be highest for the lowest CD4 category because 
%   MLE might suggest progressions that were not in the same rank order of CD4 categories

%% Find optimal parameter set
% loop over multiple start points, then compare local optima to find global
% optimum

%numStartPts = 10;
%x = zeros(length(lb),numStartPts);

initParams = [epsA;
              epsR;
              prepOut;
              artOut;
              maleActs(:,1);
              maleActs(:,2);
              maleActs(:,3);
              femaleActs(:,1);
              femaleActs(:,2);
              femaleActs(:,3);
              beta_hrHPV_val;
              beta_lrHPV_val;
              beta_lrHPV_val;
              hpv_hivMult;
              1.0;
              1.0;
              1.0;
              1.0;
              1.0;
              1.0;
              1.0;
              1.0;
              hpv_hivClear;
              hpv_hivClear;
              c3c2Mults;
              c2c1Mults;
              ones(age,1);
              kRL;
              kDR;
              kCCDet;
              kCD4(1,:,1)';
              kCD4(1,:,2)';
              kCD4(1,:,3)';
              kCD4(1,:,4)';
              kCD4(2,:,1)';
              kCD4(2,:,2)';
              kCD4(2,:,3)';
              kCD4(2,:,4)';
              kVl(1,:,1)';
              kVl(1,:,2)';
              kVl(1,:,3)';
              kVl(1,:,4)';
              kVl(2,:,1)';
              kVl(2,:,2)';
              kVl(2,:,3)';
              kVl(2,:,4)';
              [0.4 ; 0.4];
              [0.5 ; 0.5];
              1.0];

%for i = 1 : numStartPts  
options = psoptimset('UseParallel' , true , 'Cache' , 'on' ,...
    'CacheTol' , 0.1 , 'CompletePoll' , 'on' , 'TolMesh' , 0.1, ...
    'Display','iter','PlotFcn',@psplotbestf);
%x(:,i) = patternsearch(@calibratorAll, initParams , [] , [] , [] , [] , lb , ub , [] , options);
x = patternsearch(@calibratorAll, initParams , [] , [] , [] , [] , lb , ub , [] , options);
%disp('Pattern search iteration' , i)
profile viewer
    
    %initParams = lb + rand(size(lb)).*(ub - lb);
%end

%% Save calibrated parameters
file = 'HPV_calib_28Jan19.dat';
paramDir = [pwd , '\Params\'];
csvwrite([paramDir, file] , x)
