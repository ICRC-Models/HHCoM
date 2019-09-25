%% Load and save parameters, calibration data, indices, and matrices
clear all; clc

tic

paramDir = [pwd , '\Params\'];

%% Load parameters previously saved in "calibratedParams" file
load([paramDir , 'calibratedParams'] , 'popInit' , 'riskDistM' , 'riskDistF' , ...
    'mue' , 'fertility' , 'fertility2' , 'partnersM' , 'partnersF' , 'muHIV' , ...
    'maleActs' , 'femaleActs' , 'lambdaMultImm' , 'kVl' , 'circ' , 'yr' , ...
    'betaHIV_F2M' , 'betaHIV_M2F');
muHIV(11 , 2) = 0.02; % fix typo

%% Set and save general parameters
stepsPerYear = 6;
timeStep = 1 / stepsPerYear;

disease = 8;
viral = 6;
hpvVaxStates = 7;
hpvNonVaxStates = 7;
endpoints = 4;
intervens = 4;
gender = 2;
age = 80;
risk = 3;

dim = [disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , intervens , gender , age , risk];

% index retrieval function
k = cumprod([disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , intervens , gender , age]);

toInd = @(x) (x(: , 9) - 1) * k(8) + (x(: , 8) - 1) * k(7) + (x(: , 7) - 1) * k(6) + ...
    (x(: , 6) - 1) * k(5) + (x(: , 5) - 1) * k(4) + (x(: , 4) - 1) * k(3) + ...
    (x(: , 3) - 1) * k(2) + (x(: , 2) - 1) * k(1) + x(: , 1);

sumall = @(x) sum(x(:));

save(fullfile(paramDir ,'generalParams'), 'stepsPerYear' , 'timeStep' , ...
    'disease' , 'viral' , 'hpvVaxStates' , 'hpvNonVaxStates' , 'endpoints' , ...
    'intervens' , 'gender' , 'age' , 'risk' , 'dim' , 'k' , 'toInd' , ...
    'sumall');

%% Import CIN transition data from Excel
kCin1_Inf = zeros(age/5,2);
kCin2_Cin1 = zeros(age/5,2);
kCin3_Cin2 = zeros(age/5,2);
kCC_Cin3 = zeros(age/5,2);
rNormal_Inf = zeros(age/5,2);
kInf_Cin1 = zeros(age/5,2);
kCin1_Cin2 = zeros(age/5,2);
kCin2_Cin3 = zeros(age/5,2);
file = [pwd , '\Config\HPV_parameters.xlsx'];
% HPV 16/18
kCin1_Inf_orig(: , 1) = xlsread(file , 'CIN Transition' , 'B5 : B20'); % HPV to CIN1
kCin2_Cin1_orig(: , 1) = xlsread(file , 'CIN Transition' , 'C5 : C20'); % CIN1 to CIN2
kCin3_Cin2_orig(: , 1) = xlsread(file , 'CIN Transition', 'D5 : D20'); %CIN2 to CIN3
kCC_Cin3_orig(: , 1) = xlsread(file , 'CIN Transition' , 'E5 : E20'); % CIN3 to unlocalized
rNormal_Inf_orig(: , 1) = xlsread(file , 'CIN Transition' , 'F5 : F20'); % HPV to Well (natural immunity)
kInf_Cin1_orig(: , 1) = xlsread(file , 'CIN Transition' , 'G5 : G20'); % CIN1 to HPV
kCin1_Cin2_orig(: , 1) = xlsread(file , 'CIN Transition' , 'H5 : H20'); % CIN2 to CIN1
kCin2_Cin3_orig(: , 1) = xlsread(file , 'CIN Transition' , 'I5 : I20'); % CIN3 to CIN2
% 9v type ohr
kCin1_Inf_orig(: , 2) = xlsread(file , 'CIN Transition' , 'B25 : B40');
kCin2_Cin1_orig(: , 2) = xlsread(file , 'CIN Transition' , 'C25 : C40');
kCin3_Cin2_orig(: , 2) = xlsread(file , 'CIN Transition' , 'D25 : D40');
kCC_Cin3_orig(: , 2) = xlsread(file , 'CIN Transition' , 'E25 : E40');
rNormal_Inf_orig(: , 2) = xlsread(file , 'CIN Transition' , 'F25 : F40');
kInf_Cin1_orig(: , 2) = xlsread(file , 'CIN Transition' , 'G25 : G40');
kCin1_Cin2_orig(: , 2) = xlsread(file , 'CIN Transition' , 'H25 : H40');
kCin2_Cin3_orig(: , 2) = xlsread(file , 'CIN Transition' , 'I25 : I40');
% Non-vaccine type ohr
kCin1_Inf_orig(: , 3) = xlsread(file , 'CIN Transition' , 'B45 : B60');
kCin2_Cin1_orig(: , 3) = xlsread(file , 'CIN Transition' , 'C45 : C60');
kCin3_Cin2_orig(: , 3) = xlsread(file , 'CIN Transition' , 'D45 : D60');
kCC_Cin3_orig(: , 3) = xlsread(file , 'CIN Transition' , 'E45 : E60');
rNormal_Inf_orig(: , 3) = xlsread(file , 'CIN Transition' , 'F45 : F60');
kInf_Cin1_orig(: , 3) = xlsread(file , 'CIN Transition' , 'G45 : G60');
kCin1_Cin2_orig(: , 3) = xlsread(file , 'CIN Transition' , 'H45 : H60');
kCin2_Cin3_orig(: , 3) = xlsread(file , 'CIN Transition' , 'I45 : I60');
% Weight HPV transitions and HPV incidence multiplier according to type distribution
distWeight = [0.7/(0.7+0.2) , 0.2/(0.7+0.2)];
kCin1_Inf(:,1) = sum(bsxfun(@times , kCin1_Inf_orig(:,1:2) , distWeight) , 2);
kCin2_Cin1(:,1) = sum(bsxfun(@times , kCin2_Cin1_orig(:,1:2) , distWeight) , 2);
kCin3_Cin2(:,1) = sum(bsxfun(@times , kCin3_Cin2_orig(:,1:2) , distWeight) , 2);
kCC_Cin3(:,1) = sum(bsxfun(@times , kCC_Cin3_orig(:,1:2) , distWeight) , 2);
rNormal_Inf(:,1) = sum(bsxfun(@times , rNormal_Inf_orig(:,1:2) , distWeight) , 2);
kInf_Cin1(:,1) = sum(bsxfun(@times , kInf_Cin1_orig(:,1:2) , distWeight) , 2);
kCin1_Cin2(:,1) = sum(bsxfun(@times , kCin1_Cin2_orig(:,1:2) , distWeight) , 2);
kCin2_Cin3(:,1) = sum(bsxfun(@times , kCin2_Cin3_orig(:,1:2) , distWeight) , 2);
kCin1_Inf(:,2) = kCin1_Inf_orig(:,3);
kCin2_Cin1(:,2) = kCin2_Cin1_orig(:,3);
kCin3_Cin2(:,2) = kCin3_Cin2_orig(:,3);
kCC_Cin3(:,2) = kCC_Cin3_orig(:,3);
rNormal_Inf(:,2) = rNormal_Inf_orig(:,3);
kInf_Cin1(:,2) = kInf_Cin1_orig(:,3);
kCin1_Cin2(:,2) = kCin1_Cin2_orig(:,3);
kCin2_Cin3(:,2) = kCin2_Cin3_orig(:,3);

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

%% Save demographic and behavioral parameters
mInit = popInit(: , 1); % initial male population size by age
fInit = popInit(: , 2); % initial female population size by age

riskDist = zeros(age , risk , 2);
riskDist(: , : , 1) = riskDistM;
riskDist(: , : , 2) = riskDistF;

epsA = [0.3 ; 0.3 ; 0.3];
epsR = [0.3 ; 0.3 ; 0.3];
epsA_vec = cell(size(yr , 1) - 1, 1); % save data over time interval in a cell array
epsR_vec = cell(size(yr , 1) - 1, 1);
for i = 1 : size(yr , 1) - 1          % interpolate epsA/epsR values at steps within period
    period = [yr(i) , yr(i + 1)];
    epsA_vec{i} = interp1(period , epsA(i : i + 1 , 1) , ...
        yr(i) : timeStep : yr(i + 1));
    epsR_vec{i} = interp1(period , epsR(i : i + 1 , 1) , ...
        yr(i) : timeStep : yr(i + 1));
end

save(fullfile(paramDir ,'demoBehavParams'), 'mInit' , 'fInit' , 'partnersM' , 'partnersF' , ...
    'maleActs' , 'femaleActs' , 'riskDist' , 'fertility' , 'fertility2' , 'mue');

%% Save HIV natural history parameters
betaHIVF2M = zeros(age , risk , viral);
betaHIVM2F = betaHIVF2M;
for a = 1 : age % calculate per-partnership probability of HIV transmission
    betaHIVF2M(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_F2M , maleActs(a , :)')); % HIV(-) males
    betaHIVM2F(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_M2F , femaleActs(a , :)')); % HIV(-) females
end
betaHIVM2F = permute(betaHIVM2F , [2 1 3]); % risk, age, vl
betaHIVF2M = permute(betaHIVF2M , [2 1 3]); % risk, age, vl

save(fullfile(paramDir ,'hivParams'), 'betaHIVM2F' , 'betaHIVF2M' , 'muHIV' , 'kVl');

%% Save HPV natural history parameters
perPartnerHpv = 0.0045;
perPartnerHpv_nonV = 0.5 * perPartnerHpv;

% IMMUNITY
fImm(1 : age) = 1; % all infected individuals who clear HPV get natural immunity

save(fullfile(paramDir ,'hpvParams'), 'perPartnerHpv' , 'perPartnerHpv_nonV' , 'fImm' , ...
    'kCin1_Inf' , 'kCin2_Cin1' , 'kCin3_Cin2' , 'kCC_Cin3' , ...
    'rNormal_Inf' , 'kInf_Cin1' , 'kCin1_Cin2' , 'kCin2_Cin3' , 'lambdaMultImm');

%% Save intervention parameters
condUse = 0.20; %0.5 * 0.5;
OMEGA = zeros(age , 1); % hysterectomy rate

% Maximum ART coverage
maxRateM_vec = [0.40 , 0.40]; 
maxRateF_vec = [0.55 , 0.55];
maxRateM1 = maxRateM_vec(1);
maxRateM2 = maxRateM_vec(2);
maxRateF1 = maxRateF_vec(1);
maxRateF2 = maxRateF_vec(2);

% Intervention start years
hivStartYear = 1980;
circStartYear = 1990;
vaxStartYear = 2014;

cytoSens = [0.0 , 0.57 , 0.57];
% Baseline screening algorithm
baseline.screenCover = [0.0; 0.18; 0.48; 0.48; 0.48; 0.48; 0.48];
baseline.screenAge = 36;
baseline.testSens = cytoSens;
% cryoElig = [1.0 , 0.85 , 0.75 , 0.10 , 0.10 , 0.10];
baseline.colpoRetain = 0.72;
baseline.cinTreatEff = [0.905 , 0.905 , 0.766 , 0.766 , 0.766 , 0.766 , 0.766 , 0.766]; % cryotherapy/LEEP effectiveness by HIV status
baseline.cinTreatRetain = 0.51;
baseline.cinTreatHpvPersist = 0.28; % HPV persistence with LEEP
baseline.ccTreatRetain = 0.40;

save(fullfile(paramDir ,'intervenParams'), 'circ' , 'condUse' , ...
    'maxRateM1' , 'maxRateM2' , 'maxRateF1' , 'maxRateF2' , ...
    'hivStartYear' , 'circStartYear' , 'vaxStartYear' , ...
    'baseline');

%% Import and save calibration data
file = [pwd , '\Config\Calibration_targets.xlsx'];

cinPos2014_obs = xlsread(file , 'Calibration' , 'D2 : F11'); %CIN2/CIN3 Prevalence (HIV+) 2014, by age
cinNeg2014_obs = xlsread(file , 'Calibration' , 'D12 : F21'); %CIN2/CIN3 Prevalence (HIV-) 2014, by age

hpv_hiv_obs = xlsread(file , 'Calibration' , 'D144 : F152'); % HPV Prevalence in HIV+ Women (All) 2014, by age
hpv_hivNeg_obs = xlsread(file , 'Calibration' , 'D153 : F161'); % HPV Prevalence in HIV- Women (All) 2014, by age

hpv_hivM2008_obs = xlsread(file , 'Calibration' , 'E52 : F55'); % HPV Prevalence in HIV+ Men, by age
hpv_hivMNeg2008_obs = xlsread(file , 'Calibration' , 'E56 : F59'); % HPV Prevalence in HIV- Men, by age

hivPrevM_obs = xlsread(file , 'Calibration' , 'D60 : F101'); % HIV Prevalence in Men 2003,2005,2006,2007,2008,2009, by age
hivPrevF_obs = xlsread(file , 'Calibration' , 'D102 : F143'); % HIV Prevalence in Women 2003,2005,2006,2007,2008,2009, by age
 
save(fullfile(paramDir , 'calibData'), 'cinPos2014_obs' , 'cinNeg2014_obs' , ...
    'hpv_hiv_obs' , 'hpv_hivNeg_obs' , 'hpv_hivM2008_obs' , 'hpv_hivMNeg2008_obs' , ...
    'hivPrevM_obs' , 'hivPrevF_obs')




%% Load indices
disp('Preparing indices...')
disp('This may take a while...')

%% mixInfect.m indices
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

gar = zeros(gender , age , risk , disease * viral * hpvVaxStates * hpvNonVaxStates * endpoints * intervens);
for g = 1 : gender
    for a = 1 : age
        for r = 1 : risk
            gar(g , a , r , :) = sort(toInd(allcomb(1 : disease , 1 : viral ,...
                1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , g , a , r)));
        end
    end
end

hivSus = zeros(2 , gender , age , risk , hpvVaxStates * hpvNonVaxStates * endpoints * intervens);
for d = 1 : 2
    for g = 1 : gender
        for a = 1 : age
            for r = 1 : risk
                hivSus(d , g , a , r , :) =...
                    sort(toInd(allcomb(d , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , g , a , r)));
            end
        end
    end
end

hpvVaxSus = zeros(disease, gender , age , risk , viral*hpvNonVaxStates*2);
hpvVaxImm = hpvVaxSus;
hpvVaxSusVaxd = hpvVaxSus;
hpvVaxImmVaxd = hpvVaxSus;
hpvNonVaxSus = zeros(disease , gender , age , risk , viral*hpvVaxStates*2);
hpvNonVaxImm = hpvNonVaxSus;
hpvNonVaxSusVaxd = hpvNonVaxSus;
hpvNonVaxImmVaxd = hpvNonVaxSus;
for d = 1 : disease
    for g = 1 : gender
        for a = 1 : age
            for r = 1 : risk
                hpvVaxSus(d , g , a , r , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 1 , 1 : hpvNonVaxStates , 1 , [1,3] , g , a , r)));
                hpvVaxImm(d , g , a , r , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 7 , 1 : hpvNonVaxStates , 1 , [1,3] , g , a , r)));
                hpvNonVaxSus(d , g , a , r , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 , 1 , [1,3] , g , a , r)));
                hpvNonVaxImm(d , g , a , r , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 1 : hpvNonVaxStates , 7 , 1 , [1,3] , g , a , r)));

                hpvVaxSusVaxd(d , g , a , r , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 1 , 1 : hpvNonVaxStates , 1 , [2,4] , g , a , r)));
                hpvVaxImmVaxd(d , g , a , r , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 7 , 1 : hpvNonVaxStates , 1 , [2,4] , g , a , r)));
                hpvNonVaxSusVaxd(d , g , a , r , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 , 1 , [2,4] , g , a , r)));
                hpvNonVaxImmVaxd(d , g , a , r , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 7 , 1 , [2,4] , g , a , r)));
            end
        end
    end
end

toHiv = zeros(gender , age , risk , hpvVaxStates * hpvNonVaxStates * endpoints * intervens);
for g = 1 : gender
    for a = 1 : age
        for r = 1 : risk
            toHiv(g , a , r , :) = ...
                sort(toInd(allcomb(3 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , g , a , r)));
        end
    end
end

hpvVaxInfs = zeros(gender , age , risk , disease * viral * 5 * hpvNonVaxStates * endpoints * intervens);
hpvNonVaxInfs = zeros(gender , age , risk , disease , viral * hpvVaxStates * 5 * endpoints * intervens);
for g = 1 : gender
    for a = 1 : age
        for r = 1 : risk
            hpvVaxInfs(g , a , r , :) = ...
                sort(toInd(allcomb(1 : disease , 1 : viral , 2 : 6 , ...
                1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , g , a , r)));
            hpvNonVaxInfs(g , a , r , :) = ...
                sort(toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , ...
                2 : 6 , 1 : endpoints , 1 : intervens , g , a , r)));
        end
    end
end

save([paramDir , 'mixInfectIndices'] , 'mCurr' , 'fCurr' , 'mCurrArt' , 'fCurrArt' , ...
    'gar' , 'hivSus' , 'hpvVaxSus' , 'hpvVaxImm' , 'hpvVaxSusVaxd' , 'hpvVaxImmVaxd' , ...
    'hpvNonVaxSus' , 'hpvNonVaxImm' , 'hpvNonVaxSusVaxd' , 'hpvNonVaxImmVaxd' , ...
    'toHiv' , 'hpvVaxInfs' , 'hpvNonVaxInfs')
disp('mixInfect indices loaded')

%% hiv2a.m indices
hivInds = zeros(disease , viral , gender , age , risk , hpvVaxStates * hpvNonVaxStates * endpoints * intervens);
for d = 1 : disease
    for v = 1 : viral
        for g = 1 : gender
            for a = 1 : age
                for r = 1 :risk
                    hivInds(d , v , g , a , r , :) = ...
                        sort(toInd(allcomb(d , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                        1 : endpoints , 1 : intervens , g , a , r)));
                end
            end
        end
    end
end

save([paramDir , 'hivIndices'] , 'hivInds')
disp('hiv2a indices loaded')

%% hpvCCdet.m indices
disp('Preparing indices for HPV modules...')

ccLocInds = zeros(disease , age , intervens , viral * risk);
ccRegInds = ccLocInds;
ccDistInds = ccLocInds;
cin1hpvVaxInds = zeros(disease , age , intervens , viral * hpvNonVaxStates * risk);
cin2hpvVaxInds = cin1hpvVaxInds;
cin3hpvVaxInds = cin1hpvVaxInds;
cin1hpvNonVaxInds = zeros(disease , age , intervens , viral * hpvVaxStates * risk);
cin2hpvNonVaxInds = cin1hpvNonVaxInds;
cin3hpvNonVaxInds = cin1hpvNonVaxInds;
for d = 1 : disease
    for a = 1 : age
        for p = 1 : intervens
            ccLocInds(d , a , p , :) = sort(toInd(allcomb(d , 1 : viral , 6 , 6 , 1 , p , 2 , a , 1 : risk)));
            ccRegInds(d , a , p , :) = sort(toInd(allcomb(d , 1 : viral , 6 , 6 , 2 , p , 2 , a , 1 : risk)));
            ccDistInds(d , a , p , :) = sort(toInd(allcomb(d , 1 : viral , 6 , 6 , 3 , p , 2 , a , 1 : risk)));

            cin1hpvVaxInds(d , a , p , :) = sort(toInd(allcomb(d , 1 : viral , 3 , 1 : hpvNonVaxStates , 1 , p , 2 , a , 1 : risk)));
            cin2hpvVaxInds(d , a , p , :) = sort(toInd(allcomb(d , 1 : viral , 4 , 1 : hpvNonVaxStates , 1 , p , 2 , a , 1 : risk)));
            cin3hpvVaxInds(d , a , p , :) = sort(toInd(allcomb(d , 1 : viral , 5 , 1 : hpvNonVaxStates , 1 , p , 2 , a , 1 : risk)));
            cin1hpvNonVaxInds(d , a , p , :) = sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 3 , 1 , p , 2 , a , 1 : risk)));
            cin2hpvNonVaxInds(d , a , p , :) = sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 4 , 1 , p , 2 , a , 1 : risk)));
            cin3hpvNonVaxInds(d , a , p , :) = sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 5 , 1 , p , 2 , a , 1 : risk)));
        end
    end
end

normalhpvVaxInds = zeros(disease , gender , age , intervens ,  viral * hpvNonVaxStates * risk);
immunehpvVaxInds = normalhpvVaxInds;
infhpvVaxInds = normalhpvVaxInds;
normalhpvNonVaxInds = zeros(disease , gender , age , intervens ,  viral * hpvVaxStates * risk);
immunehpvNonVaxInds = normalhpvNonVaxInds;
infhpvNonVaxInds = normalhpvNonVaxInds;
for g = 1 : gender
    for d = 1 : disease
        for a = 1 : age
            for p = 1 : intervens
                normalhpvVaxInds(d , g , a , p , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 1 , 1 : hpvNonVaxStates , 1 , p , g , a , 1 : risk)));
                normalhpvNonVaxInds(d , g , a , p , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 , 1 , p , g , a , 1 : risk)));
                
                immunehpvVaxInds(d , g , a , p , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 7 , 1 : hpvNonVaxStates , 1 , p , g , a , 1 : risk)));
                immunehpvNonVaxInds(d , g , a , p , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 7 , 1 , p , g , a , 1 : risk)));
                
                infhpvVaxInds(d , g , a , p  , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 2 , 1 : hpvNonVaxStates , 1 , p , g , a , 1 : risk)));
                infhpvNonVaxInds(d , g , a , p  , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 2 , 1 , p , g , a , 1 : risk)));
            end
        end
    end      
end

% NOT UPDATED!!!!!!!!!!!!!!!!!
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

save([paramDir , 'hpvIndices'] , 'ccLocInds' , 'ccRegInds' , 'ccDistInds' , ...
    'cin1hpvVaxInds' , 'cin2hpvVaxInds' , 'cin3hpvVaxInds' , 'cin1hpvNonVaxInds' , ...
    'cin2hpvNonVaxInds' , 'cin3hpvNonVaxInds' , 'normalhpvVaxInds' , 'immunehpvVaxInds' , ...
    'infhpvVaxInds' , 'normalhpvNonVaxInds' , 'immunehpvNonVaxInds' , 'infhpvNonVaxInds')
disp('hpvCCdet indices loaded')

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

save([paramDir , 'ageRiskInds'] , 'ageInd' , 'riskInd')




%% Make matrices
pop = spalloc(prod(dim) , 1 , prod(dim));

%% Viral load progression (by CD4 count)
disp('Building viral load progression matrix')

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

save(fullfile(paramDir ,'vlAdvancer') , 'vlAdvancer' , '-v7.3')
disp('Finished building viral load progression matrix')

%% Fertility prior to 1990

% birth indices
negMaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
negFemaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 2 , 1 , 1));
posMaleBirth = toInd(allcomb(3 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
posFemaleBirth = toInd(allcomb(3 , 1 , 1 , 1 , 1 , 1 , 2 , 1 , 1));

% fertility matrix for uninfected mothers
disp('Building fertility matrix for uninfected mothers')
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
disp('Building fertility matrix for HIV-infected mothers')
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
            valsPos = [vals; ones((length(hivInfected)*2),1) .* ( 0.5*fertility(a,d-1) )];   
            xIndsNeg = [xInds; ones(length(hivInfected),1).*negMaleBirth; ones(length(hivInfected),1).*negFemaleBirth];
            yIndsNeg = [yInds; hivInfected; hivInfected];
            valsNeg = [vals; ones((length(hivInfected)*2),1) .* ( 0.5*fertility(a,d-1) )];   
        end
    end
end
hivFertPosBirth = sparse(xIndsPos , yIndsPos , valsPos , numel(pop) , numel(pop));
hivFertNegBirth = sparse(xIndsNeg , yIndsNeg , valsNeg , numel(pop) , numel(pop));

save(fullfile(paramDir ,'fertMat') , 'fertMat')
save(fullfile(paramDir ,'hivFertMats') , 'hivFertPosBirth' , 'hivFertNegBirth')

%% Fertility from 2015 onwards

% birth indices
negMaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
negFemaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 2 , 1 , 1));
posMaleBirth = toInd(allcomb(3 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
posFemaleBirth = toInd(allcomb(3 , 1 , 1 , 1 , 1 , 1 , 2 , 1 , 1));

% fertility matrix for uninfected mothers
disp('Building fertility2 matrix for uninfected mothers')
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
disp('Building fertility2 matrix for HIV-infected mothers')
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

save(fullfile(paramDir ,'fertMat2') , 'fertMat2')
save(fullfile(paramDir ,'hivFertMats2') , 'hivFertPosBirth2' , 'hivFertNegBirth2')

%% Background deaths
disp('Building death matrix')

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

save(fullfile(paramDir ,'deathMat') , 'deathMat')
disp('Death matrix complete')

%% Make circumcision matrix before current year (use when circumcision begins in model)
disp('Building circumcision matrix')

negMaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
negCircMaleBirth = toInd(allcomb(2 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
xInds = [negCircMaleBirth; negMaleBirth];
yInds = [negMaleBirth; negMaleBirth];
vals = [ones(length(negCircMaleBirth),1) .* ( circ(1) ); ones(length(negMaleBirth),1) .* ( -circ(1) )];
circMat = sparse(xInds , yInds , vals , numel(pop) , numel(pop));

save(fullfile(paramDir ,'circMat') , 'circMat')
disp('Circumcision matrix complete')

%% Make circumcision matrix after 2030
disp('Building circumcision matrix after 2030')

negMaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
negCircMaleBirth = toInd(allcomb(2 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
xInds = [negCircMaleBirth; negMaleBirth];
yInds = [negMaleBirth; negMaleBirth];
vals = [ones(length(negCircMaleBirth),1) .* ( 0.70 ); ones(length(negMaleBirth),1) .* ( -0.70 )];
circMat2 = sparse(xInds , yInds , vals , numel(pop) , numel(pop));

save(fullfile(paramDir ,'circMat2') , 'circMat2')
disp('Circumcision matrix after 2030 complete')

toc
