% Load and save parameters, calibration data, indices, and matrices
function[stepsPerYear , timeStep , startYear , currYear , endYear , ...
    years , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , ...
    intervens , gender , age , risk , hpvTypeGroups , dim , k , toInd , ...
    sumall , annlz , ...
    ageSexDebut , mInit , fInit , partnersM , partnersF , maleActs , ...
    femaleActs , riskDist , fertility , fertility2 , mue , epsA_vec , epsR_vec , ...
    yr , ...
    hivOn , betaHIVM2F , betaHIVF2M , muHIV , kVl , kCD4 , ...
    hpvOn , perPartnerHpv_vax , perPartnerHpv_nonV , fImm , rImmune , ...
    kCin1_Inf , kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , rNormal_Inf , kInf_Cin1 , ...
    kCin1_Cin2 , kCin2_Cin3 , lambdaMultImm , hpv_hivClear , rImmuneHiv , ...
    c3c2Mults , c2c1Mults , muCC , kRL , kDR , artHpvMult , ...
    hpv_hivMult , ...
    circ , condUse , screenYrs , hpvScreenStartYear , waning , ...
    maxRateM1 , maxRateF1 , maxRateM2 , maxRateF2 , hivStartYear , ...
    circStartYear , vaxStartYear , baseline , cisnet , who , whob , ...
    circProtect , condProtect , MTCTRate , hyst , ...
    OMEGA , ...
    ccInc2011_dObs , cc_dist_dObs , cin3_dist_dObs , ...
    cin1_dist_dObs , hpv_dist_dObs , cinPos2002_dObs , cinNeg2002_dObs , ...
    hpv_hiv_dObs , hpv_hivNeg_dObs , hpv_hivM2008_dObs , hpv_hivMNeg2008_dObs , ...
    hivPrevM_dObs , hivPrevF_dObs , ...
    mCurr , fCurr , mCurrArt , fCurrArt , ...
    gar , hivSus , hpvVaxSus , hpvVaxImm , hpvNonVaxSus , hpvNonVaxImm , ...
    toHiv , vaxInds , nonVInds , hpvVaxInf , hpvNonVaxInf , ...
    hivInds , ...
    cin3hpvVaxIndsFrom , ccLochpvVaxIndsTo , ccLochpvVaxIndsFrom , ...
    ccReghpvVaxInds , ccDisthpvVaxInds , cin3hpvNonVaxIndsFrom , ...
    ccLochpvNonVaxIndsTo , ccLochpvNonVaxIndsFrom , ccReghpvNonVaxInds , ...
    ccDisthpvNonVaxInds , cin1hpvVaxInds , cin2hpvVaxInds , cin3hpvVaxInds , ...
    cin1hpvNonVaxInds , cin2hpvNonVaxInds , cin3hpvNonVaxInds , normalhpvVaxInds , ...
    immunehpvVaxInds , infhpvVaxInds , normalhpvNonVaxInds , immunehpvNonVaxInds , ...
    infhpvNonVaxInds , ...
    ageInd , riskInd , ...
    vlAdvancer , ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , deathMat , circMat , circMat2] = loadUp2(fivYrAgeGrpsOn , calibBool , pIdx , paramsSub , paramSet)

% tic

paramDir = [pwd , '/Params/'];

%% Set and save general parameters

% Time
stepsPerYear = 6;
timeStep = 1 / stepsPerYear;
startYear = 1910;
currYear = 2020;
endYear = currYear; %2015; %currYear;
years = endYear - startYear;

% Compartments
disease = 8;
viral = 6;
hpvVaxStates = 7;
hpvNonVaxStates = 7;
endpoints = 4;
intervens = 4;
gender = 2;
age = 80 / max(1,fivYrAgeGrpsOn*5);
risk = 3;

hpvTypeGroups = 2;

dim = [disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , intervens , gender , age , risk];

% Index retrieval function
k = cumprod([disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , intervens , gender , age]);

toInd = @(x) (x(: , 9) - 1) * k(8) + (x(: , 8) - 1) * k(7) + (x(: , 7) - 1) * k(6) + ...
    (x(: , 6) - 1) * k(5) + (x(: , 5) - 1) * k(4) + (x(: , 4) - 1) * k(3) + ...
    (x(: , 3) - 1) * k(2) + (x(: , 2) - 1) * k(1) + x(: , 1);

sumall = @(x) sum(x(:));

annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear));

% save(fullfile(paramDir ,'generalParams'), 'stepsPerYear' , 'timeStep' , ...
%     'startYear' , 'currYear' , 'endYear' , 'years' , ...
%     'disease' , 'viral' , 'hpvVaxStates' , 'hpvNonVaxStates' , 'endpoints' , ...
%     'intervens' , 'gender' , 'age' , 'risk' , 'hpvTypeGroups' , 'dim' , 'k' , 'toInd' , ...
%     'sumall' , 'annlz');

%% Load parameters previously saved in "calibratedParams" file
load([paramDir , 'calibratedParams'] , 'popInit' , 'riskDistM' , 'riskDistF' , ...
    'mue' , 'muHIV' , 'kVl' , 'kCD4' , 'circ' , 'yr' , ...
    'fertility' , 'fertility2' , 'rImmuneHiv' , 'muCC' , ...
    'kRL' , 'kDR' , 'hpv_hivMult' , 'circProtect' , 'condProtect' , 'MTCTRate');
muHIV(11 , 2) = 0.02; % fix typo

% Male partners per year by age and risk group
if calibBool && any(1 == pIdx)
    idx = find(1 == pIdx);
    partnersM = zeros(age , risk);
    partnersMmult = paramSet(paramsSub{idx}.inds(:));
    %rowL = paramsSub{idx}.length/3;
    %rl = paramsSub{idx}.inds(1:rowL);
    %rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    %rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    partnersM(1:2 , 1:risk) = ones(2 , risk) .* 0.00001;
    %partnersM(3:10, 1:risk) = [paramSet(rl).*paramSet(rm).*paramSet(rh) , paramSet(rm).*paramSet(rh) , paramSet(rh)];
    %partnersM(11:age , 1:risk) = ones(6,risk).*partnersM(10 , 1:risk);
    partnersM(3:age , 1:risk) = partnersM(3:age , 1:risk) .* partnersMmult;
    % partnersM(10:age , 3) = ones(7 , 1);
else 
    load([paramDir , 'calibratedParams'] , 'partnersM');
end

% Female partners per year by age and risk group
if calibBool && any(2 == pIdx)
    idx = find(2 == pIdx);
    partnersF = zeros(age , risk);
    partnersFmult = paramSet(paramsSub{idx}.inds(:));
    %rowL = paramsSub{idx}.length/3;
    %rl = paramsSub{idx}.inds(1:rowL);
    %rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    %rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    partnersF(1:2 , 1:risk) = ones(2 , risk) .* 0.00001;
    %partnersF(3:10 , 1:risk) = [paramSet(rl).*paramSet(rm).*paramSet(rh) , paramSet(rm).*paramSet(rh) , paramSet(rh)];
    %partnersF(11:age , 1:risk) = ones(6,risk).*partnersF(10 , 1:risk);
    partnersF(3:age , 1:risk) = partnersF(3:age , 1:risk) .* partnersFmult;
    % partnersF(10:age , 3) = ones(7 , 1);
else
    load([paramDir , 'calibratedParams'] , 'partnersF');
end    

% Male acts per partnership per year by age and risk group
if calibBool && any(8 == pIdx)
    idx = find(8 == pIdx);
    maleActs = zeros(age , risk);
    maleActsmult = paramSet(paramsSub{idx}.inds(:));
    %rowL = paramsSub{idx}.length/3;
    %rl = paramsSub{idx}.inds(1:rowL);
    %rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    %rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    maleActs(1:2 , 1:risk) = zeros(2 , risk);
    %maleActs(3:age , 1:risk) = [paramSet(rl) , paramSet(rm).*paramSet(rl) , paramSet(rh).*paramSet(rm).*paramSet(rl)];
    maleActs(3:age , 1:risk) = maleActs(3:age , 1:risk) .* maleActsmult;
else
    load([paramDir , 'calibratedParams'] , 'maleActs');
end

% Female acts per partnership per year by age and risk group
if calibBool && any(9 == pIdx)
    idx = find(9 == pIdx);
    femaleActs = zeros(age , risk);
    femaleActsmult = paramSet(paramsSub{idx}.inds(:));
    %rowL = paramsSub{idx}.length/3;
    %rl = paramsSub{idx}.inds(1:rowL);
    %rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    %rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    femaleActs(1:2 , 1:risk) = zeros(2 , risk);
    %femaleActs(3:age , 1: risk) = [paramSet(rl) , paramSet(rm).*paramSet(rl) , paramSet(rh).*paramSet(rm).*paramSet(rl)];
    femaleActs(3:age , 1: risk) = femaleActs(3:age , 1: risk) .* femaleActsmult;
else
    load([paramDir , 'calibratedParams'] , 'femaleActs');
end

% HPV clearance multiplier for HIV-positive persons
if calibBool && any(14 == pIdx)
    idx = find(14 == pIdx);
    hpv_hivClear = zeros(4 , 1);
    hpv_hivClear(1,1) = paramSet(paramsSub{idx}.inds(1));
    hpv_hivClear(2,1) = hpv_hivClear(1,1)*paramSet(paramsSub{idx}.inds(2));
    hpv_hivClear(3,1) = hpv_hivClear(2,1)*paramSet(paramsSub{idx}.inds(3));
    hpv_hivClear(4,1) = hpv_hivClear(3,1)*paramSet(paramsSub{idx}.inds(4));
else
    load([paramDir , 'calibratedParams'] , 'hpv_hivClear');
end

% CIN2 to CIN3 progression multiplier for HIV-positive women
if calibBool && any(15 == pIdx)
    idx = find(15 == pIdx);
    c3c2Mults = zeros(4 , 1);
    c3c2Mults(4,1) = paramSet(paramsSub{idx}.inds(3));
    c3c2Mults(3,1) = c3c2Mults(4,1)*paramSet(paramsSub{idx}.inds(2));
    c3c2Mults(2,1) = c3c2Mults(3,1)*paramSet(paramsSub{idx}.inds(1));
else
    load([paramDir , 'calibratedParams'] , 'c3c2Mults');
end

% CIN1 to CIN2 progression multiplier for HIV-positive women
if calibBool && any(16 == pIdx)
    idx = find(16 == pIdx);
    c2c1Mults = zeros(4 , 1);
    c2c1Mults(4,1) = paramSet(paramsSub{idx}.inds(3));
    c2c1Mults(3,1) = c3c2Mults(4,1)*paramSet(paramsSub{idx}.inds(2));
    c2c1Mults(2,1) = c3c2Mults(3,1)*paramSet(paramsSub{idx}.inds(1));
else
    load([paramDir , 'calibratedParams'] , 'c2c1Mults');
end

% Natural immunity multiplier
if calibBool && any(18 == pIdx)
    idx = find(18 == pIdx);
    lambdaMultImm = zeros(age , 1);
    lambdaMultImmmult = paramSet(paramsSub{idx}.inds(:));
    lambdaMultImm = lambdaMultImm .* lambdaMultImmmult;
else
    load([paramDir , 'calibratedParams'] , 'lambdaMultImm');
end

%% Import HPV/CIN/CC transition data from Excel
% file = [pwd , '/Config/HPV_parameters.xlsx'];
% 
% ageTrends = xlsread(file , 'CIN transitions by type' , 'C38 : J40');
% 
% kCin1_Inf_orig = zeros(1,2);
% kCin2_Cin1_orig = zeros(1,2);
% kCin3_Cin2_orig = zeros(1,2);
% kCC_Cin3_orig = zeros(1,2);
% rNormal_Inf_orig = zeros(1,2);
% kInf_Cin1_orig = zeros(1,2);
% kCin1_Cin2_orig = zeros(1,2);
% kCin2_Cin3_orig = zeros(1,2);
% % 9v HPV types
% kCin1_Inf_orig(1 , 1) = xlsread(file , 'CIN transitions by type' , 'C25'); % HPV to CIN1, ages 10-24
% kCin2_Cin1_orig(1 , 1) = xlsread(file , 'CIN transitions by type' , 'D25'); % CIN1 to CIN2
% kCin3_Cin2_orig(1, 1) = xlsread(file , 'CIN transitions by type', 'E25'); % CIN2 to CIN3
% kCC_Cin3_orig(1 , 1) = xlsread(file , 'CIN transitions by type' , 'F25'); % CIN3 to unlocalized
% rNormal_Inf_orig(1 , 1) = xlsread(file , 'CIN transitions by type' , 'G25'); % HPV to Well (natural immunity)
% kInf_Cin1_orig(1 , 1) = xlsread(file , 'CIN transitions by type' , 'H25'); % CIN1 to HPV
% kCin1_Cin2_orig(1 , 1) = xlsread(file , 'CIN transitions by type' , 'I25'); % CIN2 to CIN1
% kCin2_Cin3_orig(1 , 1) = xlsread(file , 'CIN transitions by type' , 'J25'); % CIN3 to CIN2
% % Non-9v HPV types
% kCin1_Inf_orig(1 , 2) = xlsread(file , 'CIN transitions by type' , 'M25');  % HPV to CIN1, ages 10-24
% kCin2_Cin1_orig(1 , 2) = xlsread(file , 'CIN transitions by type' , 'N25'); % CIN1 to CIN2
% kCin3_Cin2_orig(1 , 2) = xlsread(file , 'CIN transitions by type' , 'O25'); %CIN2 to CIN3
% kCC_Cin3_orig(1 , 2) = xlsread(file , 'CIN transitions by type' , 'P25'); % CIN3 to unlocalized
% rNormal_Inf_orig(1 , 2) = xlsread(file , 'CIN transitions by type' , 'Q25'); % HPV to Well (natural immunity)
% kInf_Cin1_orig(1 , 2) = xlsread(file , 'CIN transitions by type' , 'R25'); % CIN1 to HPV
% kCin1_Cin2_orig(1 , 2) = xlsread(file , 'CIN transitions by type' , 'S25'); % CIN2 to CIN1
% kCin2_Cin3_orig(1 , 2) = xlsread(file , 'CIN transitions by type' , 'T25'); % CIN3 to CIN2
% 
% save(fullfile(paramDir ,'hpvCinCCTrans'), 'kCin1_Inf_orig' , 'kCin2_Cin1_orig' , ...
%     'kCin3_Cin2_orig' , 'kCC_Cin3_orig' , 'rNormal_Inf_orig' , 'kInf_Cin1_orig' , ...
%     'kCin1_Cin2_orig' , 'kCin2_Cin3_orig' , 'ageTrends');

%% Set HPV/CIN/CC transitions by age group based on pre-saved or calibrated values

% Set transitions in first age group from pre-saved or calibrated values
load([paramDir , 'hpvCinCCTrans'] , 'kCin1_Inf_orig' , 'kCin2_Cin1_orig' , ...
    'kCin3_Cin2_orig' , 'kCC_Cin3_orig' , 'rNormal_Inf_orig' , 'kInf_Cin1_orig' , ...
    'kCin1_Cin2_orig' , 'kCin2_Cin3_orig' , 'ageTrends');

kCin1_Inf = zeros(16,2);
kCin2_Cin1 = zeros(16,2);
kCin3_Cin2 = zeros(16,2);
kCC_Cin3 = zeros(16,2);
rNormal_Inf = zeros(16,2);
kInf_Cin1 = zeros(16,2);
kCin1_Cin2 = zeros(16,2);
kCin2_Cin3 = zeros(16,2);

% HPV to CIN1, ages 10-24
if calibBool && any(27 == pIdx)
    idx = find(27 == pIdx);
    kCin1_InfMult = paramSet(paramsSub{idx}.inds(:));
    kCin1_Inf(1 : 5 , 1) = kCin1_Inf_orig(1 , 1) * kCin1_InfMult(1);
    kCin1_Inf(1 : 5 , 2) = kCin1_Inf_orig(1 , 2) * kCin1_InfMult(2);
else
    kCin1_Inf(1 : 5 , 1) = kCin1_Inf_orig(1 , 1);
    kCin1_Inf(1 : 5 , 2) = kCin1_Inf_orig(1 , 2)*1.6;
end

% CIN1 to CIN2, ages 10-24
if calibBool && any(28 == pIdx)
    idx = find(28 == pIdx);
    kCin2_Cin1Mult = paramSet(paramsSub{idx}.inds(:));
    kCin2_Cin1(1 : 5 , 1) = kCin2_Cin1_orig(1 , 1) * kCin2_Cin1Mult(1);
    kCin2_Cin1(1 : 5 , 2) = kCin2_Cin1_orig(1 , 2) * kCin2_Cin1Mult(2);
else
    kCin2_Cin1(1 : 5 , 1) = kCin2_Cin1_orig(1 , 1);
    kCin2_Cin1(1 : 5 , 2) = kCin2_Cin1_orig(1 , 2)*1.7;
end

% CIN2 to CIN3, ages 10-24
if calibBool && any(29 == pIdx)
    idx = find(29 == pIdx);
    kCin3_Cin2Mult = paramSet(paramsSub{idx}.inds(:));
    kCin3_Cin2(1 : 5 , 1) = kCin3_Cin2_orig(1 , 1) * kCin3_Cin2Mult(1);
    kCin3_Cin2(1 : 5 , 2) = kCin3_Cin2_orig(1 , 2) * kCin3_Cin2Mult(2);
else
    kCin3_Cin2(1 : 5 , 1) = kCin3_Cin2_orig(1 , 1);
    kCin3_Cin2(1 : 5 , 2) = kCin3_Cin2_orig(1 , 2)*1.7;
end

% CIN3 to unlocalized cancer, ages 10-24
if calibBool && any(30 == pIdx)
    idx = find(30 == pIdx);
    kCC_Cin3Mult = paramSet(paramsSub{idx}.inds(:));
    kCC_Cin3(1 : 5 , 1) = kCC_Cin3_orig(1 , 1) * kCC_Cin3Mult(1);
    kCC_Cin3(1 : 5 , 2) = kCC_Cin3_orig(1 , 2) * kCC_Cin3Mult(2);
else
    kCC_Cin3(1 : 5 , 1) = kCC_Cin3_orig(1 , 1);
    kCC_Cin3(1 : 5 , 2) = kCC_Cin3_orig(1 , 2)*1.8;
end

% HPV to Well (or natural immunity), ages 10-24
if calibBool && any(31 == pIdx)
    idx = find(31 == pIdx);
    rNormal_InfMult = paramSet(paramsSub{idx}.inds(:));
    rNormal_Inf(1 : 5 , 1) = rNormal_Inf_orig(1 , 1) * rNormal_InfMult(1);
    rNormal_Inf(1 : 5 , 2) = rNormal_Inf_orig(1 , 2) * rNormal_InfMult(2);
else
    rNormal_Inf(1 : 5 , 1) = rNormal_Inf_orig(1 , 1)*1.25;
    rNormal_Inf(1 : 5 , 2) = rNormal_Inf_orig(1 , 2)*0.75;
end

% CIN1 to HPV, ages 10-24
if calibBool && any(32 == pIdx)
    idx = find(32 == pIdx);
    kInf_Cin1Mult = paramSet(paramsSub{idx}.inds(:));
    kInf_Cin1(1 : 5 , 1) = kInf_Cin1_orig(1 , 1) * kInf_Cin1Mult(1);
    kInf_Cin1(1 : 5 , 2) = kInf_Cin1_orig(1 , 2) * kInf_Cin1Mult(2);
else
    kInf_Cin1(1 : 5 , 1) = kInf_Cin1_orig(1 , 1);
    kInf_Cin1(1 : 5 , 2) = kInf_Cin1_orig(1 , 2)*0.4;
end

% CIN2 to CIN1, ages 10-24
if calibBool && any(33 == pIdx)
    idx = find(33 == pIdx);
    kCin1_Cin2Mult = paramSet(paramsSub{idx}.inds(:));
    kCin1_Cin2(1 : 5 , 1) = kCin1_Cin2_orig(1 , 1) * kCin1_Cin2Mult(1);
    kCin1_Cin2(1 : 5 , 2) = kCin1_Cin2_orig(1 , 2) * kCin1_Cin2Mult(2);
else
    kCin1_Cin2(1 : 5 , 1) = kCin1_Cin2_orig(1 , 1);
    kCin1_Cin2(1 : 5 , 2) = kCin1_Cin2_orig(1 , 2)*0.3;
end

% CIN3 to CIN2, ages 10-24
if calibBool && any(34 == pIdx)
    idx = find(34 == pIdx);
    kCin2_Cin3Mult = paramSet(paramsSub{idx}.inds(:));
    kCin2_Cin3(1 : 5 , 1) = kCin2_Cin3_orig(1 , 1) * kCin2_Cin3Mult(1);
    kCin2_Cin3(1 : 5 , 2) = kCin2_Cin3_orig(1 , 2) * kCin2_Cin3Mult(2);
else
    kCin2_Cin3(1 : 5 , 1) = kCin2_Cin3_orig(1 , 1);
    kCin2_Cin3(1 : 5 , 2) = kCin2_Cin3_orig(1 , 2)*0.3;
end

% Apply age trends to 9v HPV transitions
kCin1_Inf(6 : 10 , 1) = kCin1_Inf(1 , 1) * ageTrends(1,1); % ages 25-49
kCin2_Cin1(6 : 10 , 1) = kCin2_Cin1(1 , 1) * ageTrends(1,2);
kCin3_Cin2(6 : 10 , 1) = kCin3_Cin2(1, 1) * ageTrends(1,3);
kCC_Cin3(6 : 10 , 1) = kCC_Cin3(1 , 1) * ageTrends(1,4);
rNormal_Inf(6 : 10 , 1) = rNormal_Inf(1 , 1) * ageTrends(1,5);
kInf_Cin1(6 : 10 , 1) = kInf_Cin1(1 , 1) * ageTrends(1,6);
kCin1_Cin2(6 : 10 , 1) = kCin1_Cin2(1 , 1) * ageTrends(1,7);
kCin2_Cin3(6 : 10 , 1) = kCin2_Cin3(1 , 1) * ageTrends(1,8);
kCin1_Inf(11 : 14 , 1) = kCin1_Inf(1 , 1) * ageTrends(2,1); % ages 50-69
kCin2_Cin1(11 : 14 , 1) = kCin2_Cin1(1 , 1) * ageTrends(2,2);
kCin3_Cin2(11 : 14 , 1) = kCin3_Cin2(1 , 1) * ageTrends(2,3);
kCC_Cin3(11 : 14 , 1) = kCC_Cin3(1 , 1) * ageTrends(2,4);
rNormal_Inf(11 : 14 , 1) = rNormal_Inf(1 , 1) * ageTrends(2,5);
kInf_Cin1(11 : 14 , 1) = kInf_Cin1(1 , 1) * ageTrends(2,6);
kCin1_Cin2(11 : 14 , 1) = kCin1_Cin2(1 , 1) * ageTrends(2,7);
kCin2_Cin3(11 : 14 , 1) = kCin2_Cin3(1 , 1) * ageTrends(2,8);
kCin1_Inf(15 : 16 , 1) = kCin1_Inf(1 , 1) * ageTrends(3,1); % ages 70-79
kCin2_Cin1(15 : 16 , 1) = kCin2_Cin1(1 , 1) * ageTrends(3,2);
kCin3_Cin2(15 : 16 , 1) = kCin3_Cin2(1 , 1) * ageTrends(3,3);
kCC_Cin3(15 : 16 , 1) = kCC_Cin3(1 , 1) * ageTrends(3,4);
rNormal_Inf(15 : 16 , 1) = rNormal_Inf(1 , 1) * ageTrends(3,5);
kInf_Cin1(15 : 16 , 1) = kInf_Cin1(1 , 1) * ageTrends(3,6);
kCin1_Cin2(15 : 16 , 1) = kCin1_Cin2(1 , 1) * ageTrends(3,7);
kCin2_Cin3(15 : 16 , 1) = kCin2_Cin3(1 , 1) * ageTrends(3,8);

% Apply age trends to non-9v HPV transitions
kCin1_Inf(6 : 10 , 2) = kCin1_Inf(1 , 2) * ageTrends(1,1); % ages 25-49
kCin2_Cin1(6 : 10 , 2) = kCin2_Cin1(1 , 2) * ageTrends(1,2);
kCin3_Cin2(6 : 10 , 2) = kCin3_Cin2(1, 2) * ageTrends(1,3);
kCC_Cin3(6 : 10 , 2) = kCC_Cin3(1 , 2) * ageTrends(1,4);
rNormal_Inf(6 : 10 , 2) = rNormal_Inf(1 , 2) * ageTrends(1,5);
kInf_Cin1(6 : 10 , 2) = kInf_Cin1(1 , 2) * ageTrends(1,6);
kCin1_Cin2(6 : 10 , 2) = kCin1_Cin2(1 , 2) * ageTrends(1,7);
kCin2_Cin3(6 : 10 , 2) = kCin2_Cin3(1 , 2) * ageTrends(1,8);
kCin1_Inf(11 : 14 , 2) = kCin1_Inf(1 , 2) * ageTrends(2,1); % ages 50-69
kCin2_Cin1(11 : 14 , 2) = kCin2_Cin1(1 , 2) * ageTrends(2,2);
kCin3_Cin2(11 : 14 , 2) = kCin3_Cin2(1 , 2) * ageTrends(2,3);
kCC_Cin3(11 : 14 , 2) = kCC_Cin3(1 , 2) * ageTrends(2,4);
rNormal_Inf(11 : 14 , 2) = rNormal_Inf(1 , 2) * ageTrends(2,5);
kInf_Cin1(11 : 14 , 2) = kInf_Cin1(1 , 2) * ageTrends(2,6);
kCin1_Cin2(11 : 14 , 2) = kCin1_Cin2(1 , 2) * ageTrends(2,7);
kCin2_Cin3(11 : 14 , 2) = kCin2_Cin3(1 , 2) * ageTrends(2,8);
kCin1_Inf(15 : 16 , 2) = kCin1_Inf(1 , 2) * ageTrends(3,1); % ages 70-79
kCin2_Cin1(15 : 16 , 2) = kCin2_Cin1(1 , 2) * ageTrends(3,2);
kCin3_Cin2(15 : 16 , 2) = kCin3_Cin2(1 , 2) * ageTrends(3,3);
kCC_Cin3(15 : 16 , 2) = kCC_Cin3(1 , 2) * ageTrends(3,4);
rNormal_Inf(15 : 16 , 2) = rNormal_Inf(1 , 2) * ageTrends(3,5);
kInf_Cin1(15 : 16 , 2) = kInf_Cin1(1 , 2) * ageTrends(3,6);
kCin1_Cin2(15 : 16 , 2) = kCin1_Cin2(1 , 2) * ageTrends(3,7);
kCin2_Cin3(15 : 16 , 2) = kCin2_Cin3(1 , 2) * ageTrends(3,8);

% % if calibBool && any(25 == pIdx) % CJB note: old code
% %     idx = find(25 == pIdx);
% %     kProgrsMult = paramSet(paramsSub{idx}.inds(:));
% % else
% %     kProgrsMult = 1.0;
% % end
% % 
% % if calibBool && any(26 == pIdx)
% %     idx = find(26 == pIdx);
% %     kRegrsMult = paramSet(paramsSub{idx}.inds(:));
% % else
% %     kRegrsMult = 1.0;
% % end
% % kCin1_Inf = [kCin1_Inf , kCin1_Inf .* kProgrsMult];
% % kCin2_Cin1 = [kCin2_Cin1 , kCin2_Cin1 .* kProgrsMult];
% % kCin3_Cin2 = [kCin3_Cin2 , kCin3_Cin2 .* kProgrsMult];
% % kCC_Cin3 = [kCC_Cin3 , kCC_Cin3 .* kProgrsMult];
% % rNormal_Inf = [rNormal_Inf , rNormal_Inf .* kRegrsMult];
% % kInf_Cin1 = [kInf_Cin1 , kInf_Cin1 .* kRegrsMult];
% % kCin1_Cin2 = [kCin1_Cin2 , kCin1_Cin2 .* kRegrsMult];
% % kCin2_Cin3 = [kCin2_Cin3 , kCin2_Cin3 .* kRegrsMult];

%%
if ~fivYrAgeGrpsOn
    %% Convert 5-year age groups to 1-year age groups

    % Divide popInit age groups equally into five
    popInit_orig = popInit;
    [ageDim, valDim] = size(popInit_orig);
    popInit = zeros(ageDim*5 , valDim);
    for i = 1 : ageDim
        popInit(((i-1)*5+1) : i*5 , :) = ones(5 , valDim) .* (popInit_orig(i , :)./5);
    end

    % Replicate rates across single age groups for other variables
    vars5To1_nms = {'riskDistM' , 'riskDistF' , 'mue' , 'fertility' , 'fertility2' , ...
                 'partnersM' , 'partnersF' , 'muHIV' , 'maleActs' , 'femaleActs' , 'kCin1_Inf' , ...
                 'kCin2_Cin1' , 'kCin3_Cin2' , 'kCC_Cin3' , 'rNormal_Inf' , 'kInf_Cin1' , ...
                 'kCin1_Cin2' , 'kCin2_Cin3' , 'lambdaMultImm'};
    vars5To1_vals = {riskDistM , riskDistF , mue , fertility , fertility2 , ...
                 partnersM , partnersF , muHIV , maleActs , femaleActs , kCin1_Inf , ...
                 kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , rNormal_Inf , kInf_Cin1 , ...
                 kCin1_Cin2 , kCin2_Cin3 , lambdaMultImm};    
    for j = 1 : length(vars5To1_vals)
        valsA1 = age5To1(vars5To1_vals{j});
        assignin('base', vars5To1_nms{j} , valsA1);
    end
end

%% Save demographic and behavioral parameters
ageSexDebut = (10/max(1 , fivYrAgeGrpsOn*5)+1);

% Initial population size
mInit = popInit(: , 1); % initial male population size by age
fInit = popInit(: , 2); % initial female population size by age

% Risk distribution
riskDist = zeros(age , risk , 2);
riskDist(: , : , 1) = riskDistM;
riskDist(: , : , 2) = riskDistF;

% Mixing by age group
if calibBool && any(6 == pIdx);
    idx = find(6 == pIdx);
    epsA = zeros(3 , 1);
    %epsA = paramSet(paramsSub{idx}.inds(:));
    epsA = ones(3,1).*paramSet(paramsSub{idx}.inds(:));
else
    epsA = [0.3 ; 0.3 ; 0.3];
end
% Mixing by risk group
if calibBool && any(7 == pIdx);
    idx = find(7 == pIdx);
    epsR = zeros(3 , 1);
    %epsR = paramSet(paramsSub{idx}.inds(:));
    epsR = ones(3,1).*paramSet(paramsSub{idx}.inds(:));
else
    epsR = [0.3 ; 0.3 ; 0.3];
end
epsA_vec = cell(size(yr , 1) - 1, 1); % save data over time interval in a cell array
epsR_vec = cell(size(yr , 1) - 1, 1);
for i = 1 : size(yr , 1) - 1          % interpolate epsA/epsR values at steps within period
    period = [yr(i) , yr(i + 1)];
    epsA_vec{i} = interp1(period , epsA(i : i + 1 , 1) , ...
        yr(i) : timeStep : yr(i + 1));
    epsR_vec{i} = interp1(period , epsR(i : i + 1 , 1) , ...
        yr(i) : timeStep : yr(i + 1));
end

% save(fullfile(paramDir ,'demoBehavParams'), 'ageSexDebut' , 'mInit' , 'fInit' , 'partnersM' , 'partnersF' , ...
%     'maleActs' , 'femaleActs' , 'riskDist' , 'fertility' , 'fertility2' , 'mue' , ...
%     'epsA_vec' , 'epsR_vec' , 'yr');

%% Save HIV natural history parameters
hivOn = 1; % bool to turn HIV on or off although model calibrated for HIV on

% HIV tranmission rate
analProp = [0 , 0; 0 , 0; 0 ,0]; % [risk x gender]; proportion practicing anal sex (zero)
vagTransM = 8 / 10 ^ 4 * ones(size(analProp , 1) , 1);
vagTransF = 4 / 10 ^ 4 * ones(size(analProp , 1) , 1);
transM = vagTransM .* (1 - analProp(: , 1));
transF = vagTransF .* (1 - analProp(: , 2));
betaHIV_F2M = bsxfun(@times , [7 1 5.8 6.9 11.9 0.04; 7 1 5.8 6.9 11.9 0.04; 7 1 5.8 6.9 11.9 0.04] , transF);
betaHIV_M2F = bsxfun(@times , [7 1 5.8 6.9 11.9 0.04; 7 1 5.8 6.9 11.9 0.04; 7 1 5.8 6.9 11.9 0.04] , transM);
betaHIVF2M = zeros(age , risk , viral);
betaHIVM2F = betaHIVF2M;
for a = 1 : age % calculate per-partnership probability of HIV transmission
    % force of infection: females to infect HIV-negative males, 
    % affected by betaHIV_F2M (probability of transmission from female (receptive) to male(insertive) based on female's disease state), and number of male acts
    betaHIVF2M(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_F2M , maleActs(a , :)'));
    % force of infection: males to infect HIV-negative females,
    % affected by betaHIV_M2F (probability of transmission from male (insertive) to female (receptive) based on male's disease state), and number of female acts
    betaHIVM2F(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_M2F , femaleActs(a , :)')); 
end
betaHIVM2F = permute(betaHIVM2F , [2 1 3]); % risk, age, vl
betaHIVF2M = permute(betaHIVF2M , [2 1 3]); % risk, age, vl

%save(fullfile(paramDir ,'hivParams'), 'hivOn' , 'betaHIVM2F' , 'betaHIVF2M' , 'muHIV' , 'kVl' , 'kCD4');

%% Save HPV natural history parameters
hpvOn = 1; % bool to turn HPV on or off although model not set up for HPV to be off

% HPV tranmission rates
if calibBool && any(10 == pIdx)
    idx = find(10 == pIdx);
    perPartnerHpv_vax = paramSet(paramsSub{idx}.inds(:));
else
    perPartnerHpv_vax = 0.0045;
end

if calibBool && any(11 == pIdx)
    idx = find(11 == pIdx);
    perPartnerHpv_nonV = paramSet(paramsSub{idx}.inds(:));
else
    perPartnerHpv_nonV = perPartnerHpv_vax;
end

% Immunity
rImmune = 0.024; % for HPV16, Johnson (2012)
fImm(1 : age) = 1; % all infected individuals who clear HPV get natural immunity

% HIV and ART multipliers
if calibBool && any(21 == pIdx)
    idx = find(21 == pIdx);
    artHpvMult = paramSet(paramsSub{idx}.inds(:));
else
    artHpvMult = 1.0;
end

% save(fullfile(paramDir ,'hpvParams'), 'hpvOn' , 'perPartnerHpv_vax' , 'perPartnerHpv_nonV' , ...
%     'fImm' , 'rImmune' , 'kCin1_Inf' , 'kCin2_Cin1' , 'kCin3_Cin2' , 'kCC_Cin3' , ...
%     'rNormal_Inf' , 'kInf_Cin1' , 'kCin1_Cin2' , 'kCin2_Cin3' , 'lambdaMultImm' , ...
%     'hpv_hivClear' , 'rImmuneHiv' , 'c3c2Mults' , 'c2c1Mults' , 'muCC' , ...
%     'kRL' , 'kDR' , 'artHpvMult' , 'hpv_hivMult');

%% Save intervention parameters
% Condom use
if calibBool && any(5 == pIdx);
    idx = find(5 == pIdx);
    condUse = paramSet(paramsSub{idx}.inds(:));
else
    if fivYrAgeGrpsOn
        condUse = 0.5 * 0.5;
    else
        condUse = 0.20;
    end
end

% Background hysterectomy ********NOT UPDATED!!!!!!!!!!!!!!!!!
hyst = 0; % bool to turn background hysterectomy on or off
OMEGA = zeros(age , 1); % hysterectomy rate

% Maximum ART coverage
artOutMult = 1.0; %0.95;
maxRateM1 = 0.40*artOutMult;
maxRateM2 = 0.729*artOutMult;
maxRateF1 = 0.55*artOutMult;
maxRateF2 = 0.729*artOutMult;

% Intervention start years
hivStartYear = 1980;
circStartYear = 1990;
vaxStartYear = 2014;

% Vaccination
waning = 0;    % bool to turn waning on or off

% Screening timeframe
screenYrs = [2000; 2003; 2016; currYear; 2023; 2030; 2045];
hpvScreenStartYear = screenYrs(1);

% Screening test sensitivities
cytoSens = [0.0 , 0.57 , 0.57]; % pap smear
hpvSens = [0.0 , 0.881 , 0.881]; % careHPV
hpvSensWHO = [0.0 , 0.90 , 0.94]; % HPV test

% Baseline screening algorithm
baseline.screenCover = [0.0; 0.18; 0.48; 0.48; 0.48; 0.48; 0.48];
baseline.diseaseInds = [1 : disease];
baseline.screenAge = [35/max(1 , fivYrAgeGrpsOn*5)+1];
baseline.screenAgeMults = [1.0 / max(1 , fivYrAgeGrpsOn*5)];
baseline.testSens = cytoSens;
% cryoElig = [1.0 , 0.85 , 0.75 , 0.10 , 0.10 , 0.10];
baseline.colpoRetain = 0.72;
baseline.cinTreatEff = [0.905 , 0.905 , 0.766 , 0.766 , 0.766 , 0.766 , 0.766 , 0.766]; % cryotherapy/LEEP effectiveness by HIV status
baseline.cinTreatRetain = 0.51;
baseline.cinTreatHpvPersist = 0.28; % HPV persistence with LEEP
baseline.ccTreatRetain = 0.40;
baseline.screenCover_vec = cell(size(screenYrs , 1) - 1, 1); % save data over time interval in a cell array
for i = 1 : size(screenYrs , 1) - 1          % interpolate values at steps within period
    period = [screenYrs(i) , screenYrs(i + 1)];
    baseline.screenCover_vec{i} = interp1(period , baseline.screenCover(i : i + 1 , 1) , ...
        screenYrs(i) : timeStep : screenYrs(i + 1));
end

% CISNET screening algorithm
cisnet.screenCover = [0.0; 0.18; 0.48; 0.48; 0.48; 0.70; 0.90];
cisnet.screenAge = [(35/max(1 , fivYrAgeGrpsOn*5)+1) , (45/max(1 , fivYrAgeGrpsOn*5)+1)];
cisnet.screenAgeMults = [(1.0 / max(1 , fivYrAgeGrpsOn*5)) , (1.0 / max(1 , fivYrAgeGrpsOn*5))];
cisnet.testSens = hpvSens;
cisnet.colpoRetain = 0.81*0.85; % (compliance) * (CIN2+/CC correctly identified by same-day colposcopy)
cisnet.cinTreatEff = baseline.cinTreatEff;
cisnet.cinTreatRetain = 1.0;
cisnet.cinTreatHpvPersist = 0.48; % HPV persistence with cryotherapy 
cisnet.ccTreatRetain = 1.0;
cisnet.screenCover_vec = cell(size(screenYrs , 1) - 1, 1); % save data over time interval in a cell array
for i = 1 : size(screenYrs , 1) - 1          % interpolate values at steps within period
    period = [screenYrs(i) , screenYrs(i + 1)];
    cisnet.screenCover_vec{i} = interp1(period , cisnet.screenCover(i : i + 1 , 1) , ...
        screenYrs(i) : timeStep : screenYrs(i + 1));
end

% WHO screening algorithm - version a
who.screenCover = [0.0; 0.18; 0.48; 0.48; 0.48; 0.70; 0.90]; % CJB note: removed 90% screening compliance beginning in current year
who.testSens = hpvSensWHO;
who.colpoRetain = 1.0;
who.cinTreatEff = [1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0];
who.cinTreatRetain = 0.90; % treatment compliance
who.cinTreatHpvPersist = 0.0; %100% treatment efficacy 
who.ccTreatRetain = 0.90; % treatment compliance
who.screenCover_vec = cell(size(screenYrs , 1) - 1, 1); % save data over time interval in a cell array
for i = 1 : size(screenYrs , 1) - 1          % interpolate values at steps within period
    period = [screenYrs(i) , screenYrs(i + 1)];
    who.screenCover_vec{i} = interp1(period , who.screenCover(i : i + 1 , 1) , ...
        screenYrs(i) : timeStep : screenYrs(i + 1));
end

% WHO screening algorithm - version b (to apply WHO screening parameters at different ages by HIV status)
whob.screenCover = [0.0; 0.18; 0.48; 0.48; 0.48; 0.70; 0.90]; %CJB note: removed 90% screening compliance beginning in current year
whob.screenAge = [(35/max(1 , fivYrAgeGrpsOn*5)+1) , (45/max(1 , fivYrAgeGrpsOn*5)+1)];
whob.screenAgeMults = [(1.0 / max(1 , fivYrAgeGrpsOn*5)) , (1.0 / max(1 , fivYrAgeGrpsOn*5))];
whob.testSens = hpvSensWHO;
whob.colpoRetain = 1.0;
whob.cinTreatEff = [1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0];
whob.cinTreatRetain = 0.90; % treatment compliance
whob.cinTreatHpvPersist = 0.0; %100% treatment efficacy
whob.ccTreatRetain = 0.90; % treatment compliance
whob.screenCover_vec = cell(size(screenYrs , 1) - 1, 1); % save data over time interval in a cell array
for i = 1 : size(screenYrs , 1) - 1          % interpolate values at steps within period
    period = [screenYrs(i) , screenYrs(i + 1)];
    whob.screenCover_vec{i} = interp1(period , whob.screenCover(i : i + 1 , 1) , ...
        screenYrs(i) : timeStep : screenYrs(i + 1));
end

% save(fullfile(paramDir ,'intervenParams'), 'circ' , 'condUse' , ...
%     'screenYrs' , 'hpvScreenStartYear' , 'waning' , ...
%     'maxRateM1' , 'maxRateF1' , 'maxRateM2' , 'maxRateF2' , 'hivStartYear' , ...
%     'circStartYear' , 'vaxStartYear' , 'baseline' , 'cisnet' , 'who' , 'whob' , ...
%     'circProtect' , 'condProtect' , 'MTCTRate' , 'hyst' , 'OMEGA');

%% Import and save calibration data
% file = [pwd , '/Config/Calibration_targets.xlsx'];
% 
% ccInc2011_dObs(: , 1) = xlsread(file , 'Calibration' , 'D162 : D174'); % CC Incidence Rate 2011, by age
% ccInc2011_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H162 : I174');
% 
% cc_dist_dObs(: , 1) = xlsread(file , 'Calibration' , 'D183 : D184'); % CC type distribution for 9v and non-9v HPV types
% cc_dist_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H183 : I184');
% 
% cin3_dist_dObs(: , 1) = xlsread(file , 'Calibration' , 'D189 : D190'); % CIN3 type distribution for 9v and non-9v HPV types
% cin3_dist_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H189 : I190');
% 
% cin1_dist_dObs(: , 1) = xlsread(file , 'Calibration' , 'D187 : D188'); % CIN1 type distribution for 9v and non-9v HPV types
% cin1_dist_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H187 : I188');
% 
% hpv_dist_dObs(: , 1) = xlsread(file , 'Calibration' , 'D185 : D186'); % HPV type distribution for 9v and non-9v HPV types
% hpv_dist_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H185 : I186');
% 
% cinPos2002_dObs(: , 1) = xlsread(file , 'Calibration' , 'D2 : D11'); % CIN2/CIN3 Prevalence (HIV+) 2002, by age
% cinPos2002_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H2 : I11');
% cinNeg2002_dObs(: , 1) = xlsread(file , 'Calibration' , 'D12 : D21'); % CIN2/CIN3 Prevalence (HIV-) 2002, by age
% cinNeg2002_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H12 : I21');
% 
% hpv_hiv_dObs(: , 1) = xlsread(file , 'Calibration' , 'D144 : D152'); % HPV Prevalence in HIV+ Women (All) 2014, by age
% hpv_hiv_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H144 : I152');
% hpv_hivNeg_dObs(: , 1) = xlsread(file , 'Calibration' , 'D153 : D161'); % HPV Prevalence in HIV- Women (All) 2014, by age
% hpv_hivNeg_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H153 : I161');
% 
% hpv_hivM2008_dObs(: , 1) = xlsread(file , 'Calibration' , 'D52 : D55'); % HPV Prevalence in HIV+ Men, by age
% hpv_hivM2008_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H52 : I55');
% hpv_hivMNeg2008_dObs(: , 1) = xlsread(file , 'Calibration' , 'D56 : D59'); % HPV Prevalence in HIV- Men, by age
% hpv_hivMNeg2008_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H56 : I59');
% 
% hivPrevM_dObs(: , 1) = xlsread(file , 'Calibration' , 'D60 : D101'); % HIV Prevalence in Men 2003,2005,2006,2007,2008,2009, by age
% hivPrevM_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H60 : I101');
% hivPrevM_dObs(: , 4 : 5) = xlsread(file , 'Calibration' , 'E60 : F101'); % raw data
% hivPrevF_dObs(: , 1) = xlsread(file , 'Calibration' , 'D102 : D143'); % HIV Prevalence in Women 2003,2005,2006,2007,2008,2009, by age
% hivPrevF_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H102 : I143');
% hivPrevF_dObs(: , 4 : 5) = xlsread(file , 'Calibration' , 'E102 : F143'); % raw data
%  
% save(fullfile(paramDir , 'calibData'), 'ccInc2011_dObs' , 'cc_dist_dObs' , 'cin3_dist_dObs' , ...
%     'cin1_dist_dObs' , 'hpv_dist_dObs' , 'cinPos2002_dObs' , 'cinNeg2002_dObs' , ...
%     'hpv_hiv_dObs' , 'hpv_hivNeg_dObs' , 'hpv_hivM2008_dObs' , 'hpv_hivMNeg2008_dObs' , ...
%     'hivPrevM_dObs' , 'hivPrevF_dObs')

load([paramDir , 'calibData'], 'ccInc2011_dObs' , 'cc_dist_dObs' , 'cin3_dist_dObs' , ...
    'cin1_dist_dObs' , 'hpv_dist_dObs' , 'cinPos2002_dObs' , 'cinNeg2002_dObs' , ...
    'hpv_hiv_dObs' , 'hpv_hivNeg_dObs' , 'hpv_hivM2008_dObs' , 'hpv_hivMNeg2008_dObs' , ...
    'hivPrevM_dObs' , 'hivPrevF_dObs');


%% Load indices *****************************************************************************************************************************************************************************
% disp('Preparing indices...')
% disp('This may take a while...')

%% mixInfect.m indices
gar = zeros(gender , age , risk , disease * viral * hpvVaxStates * hpvNonVaxStates * endpoints * intervens);
for g = 1 : gender
    for a = 1 : age
        for r = 1 : risk
            gar(g , a , r , :) = sort(toInd(allcomb(1 : disease , 1 : viral ,...
                1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , g , a , r)));
        end
    end
end

vaxInds = zeros(gender , age , risk , disease * viral * 5 * hpvNonVaxStates * 3 * intervens);
nonVInds = zeros(gender , age , risk , disease * viral * hpvVaxStates * 5 * 3 * intervens);
for g = 1 : gender
    for a = 1 : age
        for r = 1 : risk
            vaxInds(g , a , r , :) = ...
                sort(toInd(allcomb(1 : disease , 1 : viral , 2 : 6 , ...
                1 : hpvNonVaxStates , 1 : 3 , 1 : intervens , g , a , r)));
            nonVInds(g , a , r , :) = ...
                sort(toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , ...
                2 : 6 , 1 : 3 , 1 : intervens , g , a , r)));
        end
    end
end

hpvVaxSus = zeros(disease, gender , age , risk , intervens , viral*hpvNonVaxStates*3);
hpvVaxImm = hpvVaxSus;
hpvVaxInf = hpvVaxSus;
hpvNonVaxSus = zeros(disease , gender , age , risk , intervens , viral*hpvVaxStates*3);
hpvNonVaxImm = hpvNonVaxSus;
hpvNonVaxInf = hpvNonVaxSus;
for d = 1 : disease
    for g = 1 : gender
        for a = 1 : age
            for r = 1 : risk
                for p = 1 : intervens
                    hpvVaxSus(d , g , a , r , p , :) = ...
                        sort(toInd(allcomb(d , 1 : viral , 1 , 1 : hpvNonVaxStates , 1 : 3 , p , g , a , r)));
                    hpvVaxImm(d , g , a , r , p , :) = ...
                        sort(toInd(allcomb(d , 1 : viral , 7 , 1 : hpvNonVaxStates , 1 : 3 , p , g , a , r)));
                    hpvVaxInf(d , g , a , r , p , :) = ...
                        sort(toInd(allcomb(d , 1 : viral , 2 , 1 : hpvNonVaxStates , 1 : 3 , p , g , a , r)));

                    hpvNonVaxSus(d , g , a , r , p , :) = ...
                        sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 , 1 : 3 , p , g , a , r)));
                    hpvNonVaxImm(d , g , a , r , p , :) = ...
                        sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 7 , 1 : 3 , p , g , a , r)));
                    hpvNonVaxInf(d , g , a , r , p , :) = ...
                        sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 2 , 1 : 3 , p , g , a , r)));
                end
            end
        end
    end
end

mCurr = zeros(age , risk , viral , 5 * hpvVaxStates * hpvNonVaxStates * endpoints * intervens); % 5 HIV+ disease states
fCurr = zeros(age , risk , viral , 5 * hpvVaxStates * hpvNonVaxStates * endpoints * intervens);
mCurrArt = zeros(age , risk , hpvVaxStates * hpvNonVaxStates * endpoints * intervens); % 1 HIV+ ART disease state
fCurrArt = zeros(age , risk , hpvVaxStates * hpvNonVaxStates * endpoints * intervens);
for a = 1 : age
    for r = 1 : risk
        for v = 1 : 5
            mCurr(a , r , v , :) = toInd(allcomb(3 : 7 , v , 1 : hpvVaxStates , ...
                1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 1 , a , r));
            fCurr(a , r , v , :) = toInd(allcomb(3 : 7 , v , 1 : hpvVaxStates , ...
                1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 2 , a , r));
        end
        mCurrArt(a , r , :) = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , ...
            1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 1 , a , r));
        fCurrArt(a , r , :) = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , ...
            1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 2 , a , r));
    end
end

hivSus = zeros(2 , hpvVaxStates , hpvNonVaxStates , endpoints , gender , age , risk , intervens);
for d = 1 : 2
    for h = 1 : hpvVaxStates
        for s = 1 : hpvNonVaxStates
            for x = 1 : endpoints
                for g = 1 : gender
                    for a = 1 : age
                        for r = 1 : risk
                            hivSus(d , h , s , x , g , a , r , :) = ...
                                sort(toInd(allcomb(d , 1 , h , s , x , 1 : intervens , g , a , r)));
                        end
                    end
                end
            end
        end
    end
end

toHiv = zeros(hpvVaxStates , hpvNonVaxStates , endpoints , gender , age , risk , intervens);
for h = 1 : hpvVaxStates
    for s = 1 : hpvNonVaxStates
        for x = 1 : endpoints
            for g = 1 : gender
                for a = 1 : age
                    for r = 1 : risk
                        toHiv(h , s , x , g , a , r , :) = ...
                            sort(toInd(allcomb(3 , 1 , h , s , x , 1 : intervens , g , a , r)));
                    end
                end
            end
        end
    end
end
% save([paramDir , 'mixInfectIndices'] , 'mCurr' , 'fCurr' , 'mCurrArt' , 'fCurrArt' , ...
%     'gar' , 'hivSus' , 'hpvVaxSus' , 'hpvVaxImm' , ...
%     'hpvNonVaxSus' , 'hpvNonVaxImm' , ...
%     'toHiv' , 'vaxInds' , 'nonVInds' , 'hpvVaxInf' , 'hpvNonVaxInf');
% disp('mixInfect indices loaded')

%% hivNH.m indices
hivInds = zeros(disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , gender , age , risk , intervens);
for d = 1 : disease
    for v = 1 : viral
        for h = 1 : hpvVaxStates
            for s = 1 : hpvNonVaxStates
                for x = 1 : endpoints
                    for g = 1 : gender
                        for a = 1 : age
                            for r = 1 :risk
                                hivInds(d , v , h , s , x , g , a , r , :) = ...
                                    sort(toInd(allcomb(d , v , h , s , ...
                                    x , 1 : intervens , g , a , r)));
                            end
                        end
                    end
                end
            end
        end
    end
end

% save([paramDir , 'hivIndices'] , 'hivInds')
% disp('hiv2a indices loaded')

%% hpvCCNH.m indices
% disp('Preparing indices for HPV modules...')

cin3hpvVaxIndsFrom = zeros(disease , hpvNonVaxStates , age , viral * intervens * risk * 3);
ccLochpvVaxIndsTo = cin3hpvVaxIndsFrom;
ccLochpvVaxIndsFrom = zeros(disease , hpvNonVaxStates , age , viral * intervens * risk);
ccReghpvVaxInds = ccLochpvVaxIndsFrom;
ccDisthpvVaxInds = ccLochpvVaxIndsFrom;
cin3hpvNonVaxIndsFrom = zeros(disease , hpvVaxStates , age , viral * intervens * risk * 3);
ccLochpvNonVaxIndsTo = cin3hpvNonVaxIndsFrom;
ccLochpvNonVaxIndsFrom = zeros(disease , hpvVaxStates , age , viral * intervens * risk);
ccReghpvNonVaxInds = ccLochpvNonVaxIndsFrom;
ccDisthpvNonVaxInds = ccLochpvNonVaxIndsFrom;
cin1hpvVaxInds = zeros(disease , age , viral * hpvNonVaxStates * intervens * risk * 3);
cin2hpvVaxInds = cin1hpvVaxInds;
cin3hpvVaxInds = cin1hpvVaxInds;
cin1hpvNonVaxInds = zeros(disease , age , viral * hpvVaxStates * intervens * risk * 3);
cin2hpvNonVaxInds = cin1hpvNonVaxInds;
cin3hpvNonVaxInds = cin1hpvNonVaxInds;
for d = 1 : disease
    for a = 1 : age
        for s = 1 : hpvNonVaxStates
            cin3hpvVaxIndsFrom(d , s , a , :) = sort(toInd(allcomb(d , 1 : viral , 5 , s , 1 : 3 , 1 : intervens , 2 , a , 1 : risk)));
            ccLochpvVaxIndsTo(d , s , a , :) = sort(toInd(allcomb(d , 1 : viral , 6 , s , 1 : 3 , 1 : intervens , 2 , a , 1 : risk)));
            ccLochpvVaxIndsFrom(d , s , a , :) = sort(toInd(allcomb(d , 1 : viral , 6 , s , 1 , 1 : intervens , 2 , a , 1 : risk)));
            ccReghpvVaxInds(d , s , a , :) = sort(toInd(allcomb(d , 1 : viral , 6 , s , 2 , 1 : intervens , 2 , a , 1 : risk)));
            ccDisthpvVaxInds(d , s , a , :) = sort(toInd(allcomb(d , 1 : viral , 6 , s , 3 , 1 : intervens , 2 , a , 1 : risk)));
        end
        for h = 1 : hpvVaxStates
            cin3hpvNonVaxIndsFrom(d , h , a , :) = sort(toInd(allcomb(d , 1 : viral , h , 5 , 1 : 3 , 1 : intervens , 2 , a , 1 : risk)));
            ccLochpvNonVaxIndsTo(d , h , a , :) = sort(toInd(allcomb(d , 1 : viral , h , 6 , 1 : 3 , 1 : intervens , 2 , a , 1 : risk)));
            ccLochpvNonVaxIndsFrom(d , h , a , :) = sort(toInd(allcomb(d , 1 : viral , h , 6 , 1 , 1 : intervens , 2 , a , 1 : risk)));
            ccReghpvNonVaxInds(d , h , a , :) = sort(toInd(allcomb(d , 1 : viral , h , 6 , 2 , 1 : intervens , 2 , a , 1 : risk)));
            ccDisthpvNonVaxInds(d , h , a , :) = sort(toInd(allcomb(d , 1 : viral , h , 6 , 3 , 1 : intervens , 2 , a , 1 : risk)));
        end
        cin1hpvVaxInds(d , a , :) = sort(toInd(allcomb(d , 1 : viral , 3 , 1 : hpvNonVaxStates , 1 : 3 , 1 : intervens , 2 , a , 1 : risk)));
        cin2hpvVaxInds(d , a , :) = sort(toInd(allcomb(d , 1 : viral , 4 , 1 : hpvNonVaxStates , 1 : 3 , 1 : intervens , 2 , a , 1 : risk)));
        cin3hpvVaxInds(d , a , :) = sort(toInd(allcomb(d , 1 : viral , 5 , 1 : hpvNonVaxStates , 1 : 3 , 1 : intervens , 2 , a , 1 : risk)));
        cin1hpvNonVaxInds(d , a , :) = sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 3 , 1 : 3 , 1 : intervens , 2 , a , 1 : risk)));
        cin2hpvNonVaxInds(d , a , :) = sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 4 , 1 : 3 , 1 : intervens , 2 , a , 1 : risk)));
        cin3hpvNonVaxInds(d , a , :) = sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 5 , 1 : 3 , 1 : intervens , 2 , a , 1 : risk)));
    end
end

normalhpvVaxInds = zeros(disease , gender , age , viral * hpvNonVaxStates * intervens * risk * 3);
immunehpvVaxInds = normalhpvVaxInds;
infhpvVaxInds = normalhpvVaxInds;
normalhpvNonVaxInds = zeros(disease , gender , age , viral * hpvVaxStates * intervens * risk * 3);
immunehpvNonVaxInds = normalhpvNonVaxInds;
infhpvNonVaxInds = normalhpvNonVaxInds;
for g = 1 : gender
    for d = 1 : disease
        for a = 1 : age
                normalhpvVaxInds(d , g , a , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 1 , 1 : hpvNonVaxStates , 1 : 3 , 1 : intervens , g , a , 1 : risk)));
                normalhpvNonVaxInds(d , g , a , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 , 1 : 3 , 1 : intervens , g , a , 1 : risk)));
                
                immunehpvVaxInds(d , g , a , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 7 , 1 : hpvNonVaxStates , 1 : 3 , 1 : intervens , g , a , 1 : risk)));
                immunehpvNonVaxInds(d , g , a , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 7 , 1 : 3 , 1 : intervens , g , a , 1 : risk)));
                
                infhpvVaxInds(d , g , a , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 2 , 1 : hpvNonVaxStates , 1 : 3 , 1 : intervens , g , a , 1 : risk)));
                infhpvNonVaxInds(d , g , a , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 2 , 1 : 3 , 1 : intervens , g , a , 1 : risk)));
        end
    end      
end

% BACKGROUND HYSTERECTOMY NOT UPDATED!!!!!!!!!!!!!!!!!
% hystSusInds = zeros(disease , 9 , age , viral * hpvVaxStates * intervens * risk);
% hystInds = zeros(disease , age , viral * hpvVaxStates * intervens * risk);
% for a = 1 : age
%     for d = 1 : disease
%         for h = 1 : hpvVaxStates
%             hysthpvVaxSusInds(d , h , a , :) = toInd(allcomb(d , 1 : viral , h , ...
%                 1 : hpvNonVaxStates , 1 : intervens , 2 , a , 1 : risk));
%             hysthpvVaxInds(d , a , :) = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , ...
%                 8 , 4 , 1 : intervens , 2 , a , 1 : risk));
%         end
%         for s = 1 : hpvNonVaxStates
%             hysthpvNonVaxSusInds(d , s , a , :) = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , ...
%                 s , 1 : intervens , 2 , a , 1 : risk));
%             hysthpvNonVaxInds(d , a , :) = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , ...
%                 8 , 1 : intervens , 2 , a , 1 : risk));
%         end
%     end
% end

% save([paramDir , 'hpvIndices'] , 'cin3hpvVaxIndsFrom' , 'ccLochpvVaxIndsTo' , ...
%     'ccLochpvVaxIndsFrom' , 'ccReghpvVaxInds' , 'ccDisthpvVaxInds' , ...
%     'cin3hpvNonVaxIndsFrom' , 'ccLochpvNonVaxIndsTo' , 'ccLochpvNonVaxIndsFrom' , ...
%     'ccReghpvNonVaxInds' , 'ccDisthpvNonVaxInds' , 'cin1hpvVaxInds' , ...
%     'cin2hpvVaxInds' , 'cin3hpvVaxInds' , 'cin1hpvNonVaxInds' , ...
%     'cin2hpvNonVaxInds' , 'cin3hpvNonVaxInds' , 'normalhpvVaxInds' , 'immunehpvVaxInds' , ...
%     'infhpvVaxInds' , 'normalhpvNonVaxInds' , 'immunehpvNonVaxInds' , 'infhpvNonVaxInds')
% disp('hpvCCdet indices loaded')

%% ageRisk.m indices
ageInd = zeros(gender , age , disease * viral * hpvVaxStates * hpvNonVaxStates * endpoints * intervens * risk);
riskInd = zeros(gender , age , risk , disease * viral * hpvVaxStates * hpvNonVaxStates * endpoints * intervens);
for g = 1 : gender
    for a = 1 : age
        ageInd(g , a , :) = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , ...
            1 : hpvNonVaxStates, 1 : endpoints , 1 : intervens , g , a , 1 : risk));   
        for r = 1 : risk
            riskInd(g , a , r , :) = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , ...
                1 : hpvNonVaxStates, 1 : endpoints , 1 : intervens , g , a , r));
        end
    end
end

% save([paramDir , 'ageRiskInds'] , 'ageInd' , 'riskInd')




%% Make matrices ******************************************************************************************************************************************************************************
pop = spalloc(prod(dim) , 1 , prod(dim));

%% Viral load progression (by CD4 count)
% disp('Building viral load progression matrix')

xInds = [];
yInds = [];
vals = [];
for g = 1 : gender
    for d = 3 : 7
        % for v = 1; vlAcute
        vlAcute = toInd(allcomb(d , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , g , 1 : age , 1 : risk));
        xInds = [xInds; vlAcute];
        yInds = [yInds; vlAcute];
        vals = [vals; ones(length(vlAcute),1) .* ( -kVl(g , d-2 , 1) )];

        % for v = 2 : 4; VL < 1,000 ->  1,000 - 10,000 -> 10,000 - 50,000
        for v = 2 : 4
            vlCurr = toInd(allcomb(d , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , g , 1 : age , 1 : risk));
            vlPrev = toInd(allcomb(d , (v - 1) , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , g , 1 : age , 1 : risk));
            xInds = [xInds; vlCurr; vlCurr]; % prev -> curr
            yInds = [yInds; vlPrev; vlCurr]; % curr -> (next)
            vals = [vals; ones(length(vlCurr),1) .* ( kVl(g , d-2 , v-1) ); ones(length(vlPrev),1) .* ( -kVl(g , d-2 , v) )];
        end
        
        % for v = 5; VL > 50,000
        vl_50k = toInd(allcomb(d , 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , g , 1 : age , 1 : risk));
        vl_10kto50k = toInd(allcomb(d , 4 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , g , 1 : age , 1 : risk));
        xInds = [xInds; vl_50k];
        yInds = [yInds; vl_10kto50k];
        vals = [vals; ones(length(vl_50k),1) .* ( kVl(g , d-2 , 4) )];
    end
end
vlAdvancer = sparse(xInds , yInds , vals , numel(pop) , numel(pop));

% save(fullfile(paramDir ,'vlAdvancer') , 'vlAdvancer' , '-v7.3')
% disp('Finished building viral load progression matrix')

%% Fertility prior to 1990

% birth indices
negMaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
negFemaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 2 , 1 , 1));
posMaleBirth = toInd(allcomb(3 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
posFemaleBirth = toInd(allcomb(3 , 1 , 1 , 1 , 1 , 1 , 2 , 1 , 1));

% fertility matrix for uninfected mothers
% disp('Building fertility matrix for uninfected mothers')
xInds = [];
yInds = [];
vals = [];
for a = 1 : age
    hivUninf = toInd(allcomb(1 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : 3 , ...
        1 : intervens , 2 , a , 1 : risk));
    hivPosArt = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : 3 , ...
        1 : intervens , 2 , a , 1 : risk));
    xInds = [xInds; ones(length(hivUninf),1).*negMaleBirth; ones(length(hivUninf),1).*negFemaleBirth; ...
        ones(length(hivPosArt),1).*negMaleBirth; ones(length(hivPosArt),1).*negFemaleBirth];
    yInds = [yInds; hivUninf; hivUninf; hivPosArt; hivPosArt];
    vals = [vals; ones((length(hivUninf)*2+length(hivPosArt)*2),1) .* ( 0.5*fertility(a,1) )];
end
fertMat = sparse(xInds , yInds , vals , numel(pop) , numel(pop));

% fertility matrix for infected mothers
% disp('Building fertility matrix for HIV-infected mothers')
xIndsPos = [];
yIndsPos = [];
valsPos = [];
xIndsNeg = [];
yIndsNeg = [];
valsNeg = [];
for d = 3 : 7 % hiv infected
    for v = 1 : viral % hiv infected
        for a = 1 : age
            hivInfected = toInd(allcomb(d , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : 3 , 1 : intervens , 2 , a , 1 : risk));
            xIndsPos = [xIndsPos; ones(length(hivInfected),1).*posMaleBirth; ones(length(hivInfected),1).*posFemaleBirth];
            yIndsPos = [yIndsPos; hivInfected; hivInfected];
            valsPos = [valsPos; ones((length(hivInfected)*2),1) .* ( 0.5*fertility(a,d-1) )];   
            xIndsNeg = [xIndsNeg; ones(length(hivInfected),1).*negMaleBirth; ones(length(hivInfected),1).*negFemaleBirth];
            yIndsNeg = [yIndsNeg; hivInfected; hivInfected];
            valsNeg = [valsNeg; ones((length(hivInfected)*2),1) .* ( 0.5*fertility(a,d-1) )];   
        end
    end
end
hivFertPosBirth = sparse(xIndsPos , yIndsPos , valsPos , numel(pop) , numel(pop));
hivFertNegBirth = sparse(xIndsNeg , yIndsNeg , valsNeg , numel(pop) , numel(pop));

% save(fullfile(paramDir ,'fertMat') , 'fertMat')
% save(fullfile(paramDir ,'hivFertMats') , 'hivFertPosBirth' , 'hivFertNegBirth')

%% Fertility from 2015 onwards

% birth indices
negMaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
negFemaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 2 , 1 , 1));
posMaleBirth = toInd(allcomb(3 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
posFemaleBirth = toInd(allcomb(3 , 1 , 1 , 1 , 1 , 1 , 2 , 1 , 1));

% fertility matrix for uninfected mothers
% disp('Building fertility2 matrix for uninfected mothers')
xInds = [];
yInds = [];
vals = [];
for a = 1 : age
    hivUninf = toInd(allcomb(1 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : 3 , ...
        1 : intervens , 2 , a , 1 : risk));
    hivPosArt = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : 3 , ...
        1 : intervens , 2 , a , 1 : risk));
    xInds = [xInds; ones(length(hivUninf),1).*negMaleBirth; ones(length(hivUninf),1).*negFemaleBirth; ...
        ones(length(hivPosArt),1).*negMaleBirth; ones(length(hivPosArt),1).*negFemaleBirth];
    yInds = [yInds; hivUninf; hivUninf; hivPosArt; hivPosArt];
    vals = [vals; ones((length(hivUninf)*2+length(hivPosArt)*2),1) .* ( 0.5*fertility2(a,1) )];
end
fertMat2 = sparse(xInds , yInds , vals , numel(pop) , numel(pop));

% fertility matrix for infected mothers
% disp('Building fertility2 matrix for HIV-infected mothers')
xIndsPos = [];
yIndsPos = [];
valsPos = [];
xIndsNeg = [];
yIndsNeg = [];
valsNeg = [];
for d = 3 : 7 % hiv infected
    for v = 1 : viral % hiv infected
        for a = 1 : age
            hivInfected = toInd(allcomb(d , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : 3 , 1 : intervens , 2 , a , 1 : risk));
            xIndsPos = [xInds; ones(length(hivInfected),1).*posMaleBirth; ones(length(hivInfected),1).*posFemaleBirth];
            yIndsPos = [yInds; hivInfected; hivInfected];
            valsPos = [vals; ones((length(hivInfected)*2),1) .* ( 0.5*fertility2(a,d-1) )];   
            xIndsNeg = [xInds; ones(length(hivInfected),1).*negMaleBirth; ones(length(hivInfected),1).*negFemaleBirth];
            yIndsNeg = [yInds; hivInfected; hivInfected];
            valsNeg = [vals; ones((length(hivInfected)*2),1) .* ( 0.5*fertility2(a,d-1) )];   
        end
    end
end
hivFertPosBirth2 = sparse(xIndsPos , yIndsPos , valsPos , numel(pop) , numel(pop));
hivFertNegBirth2 = sparse(xIndsNeg , yIndsNeg , valsNeg , numel(pop) , numel(pop));

% save(fullfile(paramDir ,'fertMat2') , 'fertMat2')
% save(fullfile(paramDir ,'hivFertMats2') , 'hivFertPosBirth2' , 'hivFertNegBirth2')

%% Background deaths
% disp('Building death matrix')

xInds = [];
yInds = [];
vals = [];
for a = 1 : age
    males = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
    females = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    
    xInds = [xInds; males; females];
    yInds = [yInds; males; females];
    vals = [vals; ones(length(males),1) .* ( -mue(a,1) ); ones(length(females),1) .* ( -mue(a,2) )];
end
deathMat = sparse(xInds , yInds , vals , numel(pop) , numel(pop));

% save(fullfile(paramDir ,'deathMat') , 'deathMat')
% disp('Death matrix complete')

%% Make circumcision matrix before current year (use when circumcision begins in model)
% disp('Building circumcision matrix')

negMaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
negCircMaleBirth = toInd(allcomb(2 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
xInds = [negCircMaleBirth; negMaleBirth];
yInds = [negMaleBirth; negMaleBirth];
vals = [ones(length(negCircMaleBirth),1) .* ( circ(1) ); ones(length(negMaleBirth),1) .* ( -circ(1) )];
circMat = sparse(xInds , yInds , vals , numel(pop) , numel(pop));

% save(fullfile(paramDir ,'circMat') , 'circMat')
% disp('Circumcision matrix complete')

%% Make circumcision matrix after 2030
% disp('Building circumcision matrix after 2030')

negMaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
negCircMaleBirth = toInd(allcomb(2 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
xInds = [negCircMaleBirth; negMaleBirth];
yInds = [negMaleBirth; negMaleBirth];
vals = [ones(length(negCircMaleBirth),1) .* ( 0.70 ); ones(length(negMaleBirth),1) .* ( -0.70 )];
circMat2 = sparse(xInds , yInds , vals , numel(pop) , numel(pop));

% save(fullfile(paramDir ,'circMat2') , 'circMat2')
% disp('Circumcision matrix after 2030 complete')

% toc
