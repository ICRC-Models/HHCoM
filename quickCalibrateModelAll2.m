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
load([paramDir , 'calibratedParams'])

%% Set parameter upper and lower bounds 

% KEY
%(1): perPartnerHpv, [1x1], (0.0 to 1.0)
%(2): perPartnerHpv_lr, [1x1], (0.0 to 1.0)
%(3): perPartnerHpv_nonV, [1x1], (0.0 to 1.0)
%(4:7): hpv_hivMult, [CD4x1], (0.0 to 1.0) 
%(8): kCin1_Inf, [1x1], (0.001 to 20)
%(9): kCin2_Cin1, [1x1], (0.001 to 20)
%(10): kCin3_Cin2, [1x1], (0.001 to 20)
%(11): kCC_Cin3, [1x1], (0.001 to 20)
%(12): rNormal_Inf, [1x1], (0.001 to 20)
%(13): kInf_Cin1, [1x1], (0.001 to 20)
%(14): kCin1_Cin2, [1x1], (0.001 to 20)
%(15): kCin2_Cin3, [1x1], (0.001 to 20)
%(16:19): hpv_hivClear, [CD4x1], (0.0 to 1.0)
%(20:23): rImmuneHiv, [CD4x1], (0.0 to 1.0)
%(24:27): c3c2Mults, [CD4x1], init, (all ones to x10)
%(28:31): c2c1Mults, [CD4x1], init, (all ones to x10)
%(32): kRL, [1x1], (0.0 to 1.0)
%(33): kDR, [1x1], (0.0 to 1.0)
%(34:36): kCCDet, [3x1], (0.0 to 1.0)
%(37): artHpvMult, [1x1], (1.0 to x10)
%(38:55): muCC, [6x3], (0.0001 to 0.5)
%(56:73): muCC_det, [6x3], (0.0001 to 0.5)

%% Set parameter constraints
% max HPV acquisition multiplier for HIV+ on ART <= max acquisition multiplier for HIV+
% restrict CD4 progressions to be highest for the lowest CD4 category because 
%   MLE might suggest progressions that were not in the same rank order of CD4 categories

%% Find optimal parameter set
% loop over multiple start points, then compare local optima to find global
% optimum

%numStartPts = 10;
%x = zeros(length(lb),numStartPts);

% KEY
%(1): perPartnerHpv, [1x1], (0.0 to 1.0)
%(2): perPartnerHpv_lr, [1x1], (0.0 to 1.0)
%(3): perPartnerHpv_nonV, [1x1], (0.0 to 1.0)
%(4:7): hpv_hivMult, [CD4x1], (0.0 to 1.0) 
%(8): kCin1_Inf, [1x1], (0.001 to 20)
%(9): kCin2_Cin1, [1x1], (0.001 to 20)
%(10): kCin3_Cin2, [1x1], (0.001 to 20)
%(11): kCC_Cin3, [1x1], (0.001 to 20)
%(12): rNormal_Inf, [1x1], (0.001 to 20)
%(13): kInf_Cin1, [1x1], (0.001 to 20)
%(14): kCin1_Cin2, [1x1], (0.001 to 20)
%(15): kCin2_Cin3, [1x1], (0.001 to 20)
%(16:19): hpv_hivClear, [CD4x1], (0.0 to 1.0)
%(20:23): rImmuneHiv, [CD4x1], (0.0 to 1.0)
%(24:27): c3c2Mults, [CD4x1], init, (all ones to x10)
%(28:31): c2c1Mults, [CD4x1], init, (all ones to x10)
%(32): kRL, [1x1], (0.0 to 1.0)
%(33): kDR, [1x1], (0.0 to 1.0)
%(34:36): kCCDet, [3x1], (0.0 to 1.0)
%(37): artHpvMult, [1x1], (1.0 to x10)
%(38:55): muCC, [6x3], (0.0001 to 0.5)
%(56:73): muCC_det, [6x3], (0.0001 to 0.5)
initParams = [perPartnerHpv;
              perPartnerHpv_lr; 
              perPartnerHpv_nonV; 
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
              rImmuneHiv
              c3c2Mults
              c2c1Mults
              kRL
              kDR
              kCCDet
              artHpvMult
              muCC(:,1)
              muCC(:,2)
              muCC(:,3)
              muCC_det(:,1)
              muCC_det(:,2)
              muCC_det(:,3)];
          
lb = initParams .* 0.1;
ub = initParams .* 10;
ub(1:7) = min( ub(1:7) , 1.0);
ub(16:23) = min( ub(16:23) , 1.0);
lb(24:31) = max( lb(24:31) , 1.0);
ub(32:36) = min( ub(32:36) , 1.0);
lb(37) = max( lb(37), 1.0);
ub(38:73) = min( ub(38:73) , 1.0);

%for i = 1 : numStartPts  
options = psoptimset('UseParallel' , true , 'Cache' , 'on' ,...
    'CacheTol' , 0.1 , 'CompletePoll' , 'on' , 'TolMesh' , 0.1, ...
    'Display','iter','PlotFcn',@psplotbestf);
%x(:,i) = patternsearch(@calibratorAll, initParams , [] , [] , [] , [] , lb , ub , [] , options);
x = patternsearch(@calibratorAll2, initParams , [] , [] , [] , [] , lb , ub , [] , options);
%disp('Pattern search iteration' , i)
profile viewer
    
    %initParams = lb + rand(size(lb)).*(ub - lb);
%end

%% Save calibrated parameters
file = 'HPV_calib_11Feb19_noErr.dat';
paramDir = [pwd , '\Params\'];
csvwrite([paramDir, file] , x)
