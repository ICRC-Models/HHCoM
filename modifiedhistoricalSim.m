% Historical module
% Runs simulation over the time period and time step specified by the user.

function [] = modifiedhistoricalSim(calibBool , pIdx , paramsSub , paramSet , paramSetIdx , tstep_abc , date , sympParams_in, sympRun)
%Run from the Command Window: historicalSim(0 , [] , [] , [] , [] , 0 , '17Dec19')

%% If using pattern search algorithm, uncomment the following and change the function above to historicalSim(paramSet). 
% Note: Make sure you are calculating NEGATIVE summed log-likelihood. 
% Note: If you include incidence calibration data, need to change
% likelihood function and calibration data Excel document as modeled in likeFun.
%
% calibBool = 1;
% paramSetIdx = 1;
% tstep_abc = 0;
% date = '19Nov19';
% paramDir = [pwd , '/Params/'];
% pIdx = load([paramDir , 'pIdx_patternSrch_' , date , '_0.dat']);
% [paramsAll] = genParamStruct();
% paramsSub = cell(length(pIdx),1);
% startIdx = 1;
% for s = 1 : length(pIdx)
%     paramsSub{s}.length = paramsAll{pIdx(s)}.length;
%     paramsSub{s}.inds = (startIdx : (startIdx + paramsSub{s}.length - 1));
%     startIdx = startIdx + paramsSub{s}.length;
% end


%%
%close all; clear all; clc;
tic
% profile clear;

%%  Variables/parameters to set based on your scenario

% DIRECTORY TO SAVE RESULTS
pathModifier = ['toNow_' , date , '_stochMod_' , 'recalib_07Jul23_' , num2str(paramSetIdx)]; % ***SET ME***: name for historical run output file 
 %pathModifier = 'toNow_determMod_final_artDiscontFix';
 %pathModifier = 'toNow_determMod_popFertFix';

% AGE GROUPS
fivYrAgeGrpsOn = 1; % choose whether to use 5-year or 1-year age groups

% SCREENING
screenAlgorithm = 1; % ***SET ME***: screening algorithm to use (1 for baseline, 2 for CISNET, 3 for WHOa, 4 for WHOb)
hivPosScreen = 1; % ***SET ME***: 0 applies same screening algorithm (screenAlgorithm) for all HIV states; 1 applies screenAlgorithm to HIV+ and screenAlgorithmNeg to HIV-
screenAlgorithmNeg = 1; % ***SET ME***: If hivPosScreen=1, screening algorithm to use for HIV- persons (1 for baseline, 2 for CISNET, 3 for WHOa, 4 for WHOb) 

% VACCINATION
% CLH: I commented vaxEff out since I added the vaxEff assignment in
% loadup2. 
% vaxEff = 1.0; % actually bivalent vaccine, but to avoid adding additional compartments, we use nonavalent vaccine and then reduce coverage

%Parameters for school-based vaccination regimen  % ***SET ME***: coverage for baseline vaccination of 9-year-old girls
vaxAge = [10/max(1 , fivYrAgeGrpsOn*5)];
% vaxRate = 0.16 * vaxRateAdjust; %0.86*(0.7/0.9);    % (9 year-old coverage * bivalent vaccine efficacy adjustment)
vaxG = 2;   % indices of genders to vaccinate (1 or 2 or 1,2)

%% Vaccine scale up

vaxRateAdjust = 0.7/0.9; %bivalent/quadrivalent vaccine efficacy adjustment 

gradScaleUp = 1; % ***SET ME***: 1 if you want gradual scale up of vaccination coverage

stepsPerYear = 6; % ***SET ME***: If this changes in loadup2, you need to change it here as well
timeStep = 1 / stepsPerYear; % ***SET ME***: same here

if gradScaleUp==1
    vaxRate = [0.0; 0.16; 0.31] * vaxRateAdjust; % Coverage over time (Years: [2021; 2026])
    vaxYrs = [2019; 2020; 2023];
    vaxCover_vec = cell(size(vaxYrs , 1) - 1, 1); % save data over time interval in a cell array
    for i = 1 : size(vaxYrs , 1) - 1          % interpolate values at steps within period
        period = [vaxYrs(i) , vaxYrs(i + 1)];
        vaxCover_vec{i} = interp1(period , vaxRate(i : i + 1 , 1) , ...
            vaxYrs(i) : timeStep : vaxYrs(i + 1));
    end
    vaxRate_vec = vaxCover_vec; 
else 
    vaxRate_vec = [0.16] * vaxRateAdjust;
    vaxYrs = [2020]; % for testing 2020 orig
end 

%% Save pre-loaded parameters and pre-calculated indices and matrices
[stepsPerYear , timeStep , startYear , currYear , endYear , ...
    years , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , ...
    intervens , gender , age , risk , hpvTypeGroups , dim , k , toInd, annlz , ...
    ageSexDebut , mInit , fInit , partnersM , partnersF , partnersMmult, maleActs , ...
    femaleActs , riskDist , fertility , fertility2 , fertility3 , fertility4,...
    mue , mue2 , mue3 , mue4 , mue5, epsA_vec , epsR_vec , yr , ...
    hivOn , betaHIV_mod , hiv_hpvMult, muHIV , kCD4 , ...
    hpvOn , beta_hpvVax_mod , beta_hpvNonVax_mod , fImm , rImmune , ...
    kCin1_Inf , kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , rNormal_Inf , kInf_Cin1 , ...
    kCin1_Cin2 , kCin2_Cin3 , lambdaMultImm , hpv_hivClear , rImmuneHiv , ...
    c3c2Mults , c2c1Mults , muCC , muCC_ud , muCC_d , kRL , kDR , artHpvMult , ...
    hpv_hivMult , maleHpvClearMult , ...
    condUse , screenYrs , hpvScreenStartYear , waning , ...
    artYr , maxRateM , maxRateF , ...
    artYr_vec , artM_vec , artF_vec , minLim , maxLim , ...
    circ_aVec , vmmcYr_vec , vmmc_vec , vmmcYr , vmmcRate , ...
    hivStartYear , circStartYear ,circNatStartYear , vaxStartYear , baseline , cisnet , who , whob , ...
    circProtect , condProtect , MTCTRate , hyst , OMEGA , ...
    ccInc2012_dObs , cc_dist_dObs , cin3_dist_dObs , cin1_dist_dObs , ...
    hpv_dist_dObs , cinPos2007_dObs , cin1_2010_dObs , cin2_2010_dObs, ...
    hpv_hiv_dObs , hpv_hivNeg_dObs , hpv_all_dObs , hpv_hiv2009_dObs , ...
    hivPrevF_dObs , hivPrevM_dObs , hivPrevAll_dObs, popAgeDist_dObs , totPopSize_dObs , hivCurr , ...
    gar , hivSus , hpvVaxSus , hpvVaxImm , hpvNonVaxSus , hpvNonVaxImm , ...
    toHiv , vaxInds , nonVInds , hpvVaxInf , hpvNonVaxInf , hivInds , ...
    cin3hpvVaxIndsFrom , ccLochpvVaxIndsTo , ccLochpvVaxIndsFrom , ...
    ccReghpvVaxInds , ccDisthpvVaxInds , cin3hpvNonVaxIndsFrom , ...
    ccLochpvNonVaxIndsTo , ccLochpvNonVaxIndsFrom , ccReghpvNonVaxInds , ...
    ccDisthpvNonVaxInds , cin1hpvVaxInds , cin2hpvVaxInds , cin3hpvVaxInds , ...
    cin1hpvNonVaxInds , cin2hpvNonVaxInds , cin3hpvNonVaxInds , normalhpvVaxInds , ...
    immunehpvVaxInds , infhpvVaxInds , normalhpvNonVaxInds , immunehpvNonVaxInds , ...
    infhpvNonVaxInds , ageInd , riskInd , ...
    kSymp , hystMult , ...
    hivNegNonVMMCinds , hivNegVMMCinds , vlAdvancer , ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , fertMat3 , hivFertPosBirth3 , hivFertNegBirth3 , ...
    fertMat4 , hivFertPosBirth4 , hivFertNegBirth4 , ...
    dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ...
    dFertPos3 , dFertNeg3  , dFertMat3, d_partnersMmult, riskAdj, d_riskAdj, ...
    deathMat , deathMat2 , deathMat3 , deathMat4 , deathMat5,...
    dDeathMat , dDeathMat2 , dDeathMat3 , dDeathMat4, dMue , ...
    ccLochpvVaxIndsFrom_treat , ...
    ccReghpvVaxInds_treat , ccDisthpvVaxInds_treat , vaxEff] = loadUp2(fivYrAgeGrpsOn , calibBool , pIdx , paramsSub , paramSet , paramSetIdx);

%% Modifying the rate of symptomatic detection 

kSymp = sympParams_in(1:3);
kRL = sympParams_in(4); 
kDR = sympParams_in(5); 


%% Screening
if (screenAlgorithm == 1)
    % Baseline screening algorithm
    screenAlgs{1} = baseline;
elseif (screenAlgorithm == 2)
    % CISNET screening algorithm
    screenAlgs{1} = cisnet;
elseif (screenAlgorithm == 3)
    % WHO screening algorithm - version a
    screenAlgs{1} = who;
elseif (screenAlgorithm == 4)
    % WHO screening algorithm - version b
    screenAlgs{1} = whob;
end

if hivPosScreen
    if (screenAlgorithmNeg == 1)
        % Baseline screening algorithm
        screenAlgs{2} = baseline;
    elseif (screenAlgorithmNeg == 2)
        % CISNET screening algorithm
        screenAlgs{2} = cisnet;
    elseif (screenAlgorithmNeg == 3)
        % WHO screening algorithm - version a
        screenAlgs{2} = who;
    elseif (screenAlgorithmNeg == 4)
        % WHO screening algorithm - version b
        screenAlgs{2} = whob;
    end
    screenAlgs{2}.screenCover_vec = cell(size(screenYrs , 1) - 1, 1); % save data over time interval in a cell array
    for i = 1 : size(screenYrs , 1) - 1          % interpolate dnaTestCover values at steps within period
        period = [screenYrs(i) , screenYrs(i + 1)];
        screenAlgs{2}.screenCover_vec{i} = interp1(period , screenAlgs{2}.screenCover(i : i + 1 , 1) , ...
            screenYrs(i) : timeStep : screenYrs(i + 1));
    end
    screenAlgs{1}.diseaseInds = [3 : 8];
    screenAlgs{2}.diseaseInds = [1 : 2];
else
    screenAlgs{1}.diseaseInds = [1 : disease];
end

screenAlgs{1}.screenCover_vec = cell(size(screenYrs , 1) - 1, 1); % save data over time interval in a cell array
for i = 1 : size(screenYrs , 1) - 1          % interpolate dnaTestCover values at steps within period
    period = [screenYrs(i) , screenYrs(i + 1)];
    screenAlgs{1}.screenCover_vec{i} = interp1(period , screenAlgs{1}.screenCover(i : i + 1 , 1) , ...
        screenYrs(i) : timeStep : screenYrs(i + 1));
end

% Create screening indices
numScreenAge = length(screenAlgs{1}.screenAge);
agesComb = screenAlgs{1}.screenAge;
ageMultsComb = screenAlgs{1}.screenAgeMults;
if hivPosScreen
    numScreenAge = numScreenAge + length(screenAlgs{2}.screenAge);
    agesComb = [agesComb , screenAlgs{2}.screenAge];
    ageMultsComb = [ageMultsComb , screenAlgs{2}.screenAgeMults];
end
% difference between the KZN and Kenya models in the 3 lines below.
screenAgeAll = zeros(disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , intervens , numScreenAge , risk);
screenAgeS = zeros(disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , 2 , numScreenAge , risk);
noVaxNoScreen = zeros(disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , numScreenAge , risk);
noVaxToScreen = noVaxNoScreen;
vaxNoScreen = noVaxNoScreen;
vaxToScreen = noVaxNoScreen;
noVaxToScreenTreatImm = zeros(disease , viral , numScreenAge , risk);
vaxToScreenTreatImm = noVaxToScreenTreatImm;
noVaxToScreenTreatHpv = noVaxToScreenTreatImm;
vaxToScreenTreatHpv = noVaxToScreenTreatImm;
noVaxToScreenHyst = noVaxToScreenTreatImm;
vaxToScreenHyst = noVaxToScreenTreatImm;
noVaxToScreenTreatVaxHpv = zeros(disease , viral , hpvNonVaxStates , numScreenAge , risk);
vaxToScreenTreatVaxHpv = noVaxToScreenTreatVaxHpv;
noVaxToScreenTreatNonVaxHpv = zeros(disease , viral , hpvVaxStates , numScreenAge , risk);
vaxToScreenTreatNonVaxHpv = noVaxToScreenTreatNonVaxHpv;
noVaxScreen = zeros(disease*viral*hpvVaxStates*hpvNonVaxStates*endpoints*risk , numScreenAge);
noVaxXscreen = noVaxScreen;
vaxScreen = noVaxScreen;
vaxXscreen = noVaxScreen;
noVaxToScreenCancerNoTreat = zeros(disease , viral , 3 , numScreenAge ,risk);
noVaxToScreenCancerTreat = noVaxToScreenCancerNoTreat;
noVaxToScreenCancerNegScreen = noVaxToScreenCancerNoTreat; 
vaxToScreenCancerNegScreen = noVaxToScreenCancerNoTreat; 
vaxToScreenCancerNoTreat = noVaxToScreenCancerNoTreat;
vaxToScreenCancerTreat = noVaxToScreenCancerNoTreat;
udPop = zeros(disease, viral, hpvVaxStates, hpvNonVaxStates, 3, intervens, age, risk); 
udPopNoTreat = udPop; 
udPopTreat = udPop; 
udPopHyst = udPop;

for aS = 1 : numScreenAge
    a = agesComb(aS);
    
    for d = 1 : disease
        for v = 1 : viral
            for h = 1 : hpvVaxStates
                for s = 1 : hpvNonVaxStates
                    for x = 1 : endpoints
                        for r = 1 : risk
                            screenAgeAll(d,v,h,s,x,:,aS,r) = toInd(allcomb(d , v , h , s , x , 1 : intervens , 2 , a , r)); 
                            screenAgeS(d,v,h,s,x,:,aS,r) = toInd(allcomb(d , v , h , s , x , 3 : intervens , 2 , a , r));

                            noVaxNoScreen(d,v,h,s,x,aS,r) = sort(toInd(allcomb(d , v , h , s , x , 1 , 2 , a , r)));
                            noVaxToScreen(d,v,h,s,x,aS,r) = sort(toInd(allcomb(d , v , h , s , x , 3 , 2 , a , r)));
                            vaxNoScreen(d,v,h,s,x,aS,r) = sort(toInd(allcomb(d , v , h , s , x , 2 , 2 , a , r)));
                            vaxToScreen(d,v,h,s,x,aS,r) = sort(toInd(allcomb(d , v , h , s , x , 4 , 2 , a , r)));

                            noVaxToScreenTreatImm(d,v,aS,r) = toInd(allcomb(d , v , 7 , 7 , 1 , 3 , 2 , a , r));
                            vaxToScreenTreatImm(d,v,aS,r) = toInd(allcomb(d , v , 7 , 7 , 1 , 4 , 2 , a , r));
                            noVaxToScreenTreatHpv(d,v,aS,r) = toInd(allcomb(d , v , 2 , 2 , 1 , 3 , 2 , a , r));
                            vaxToScreenTreatHpv(d,v,aS,r) = toInd(allcomb(d , v , 2 , 2 , 1 , 4 , 2 , a , r));
                            noVaxToScreenTreatVaxHpv(d,v,s,aS,r) = toInd(allcomb(d , v , 2 , s , 1 , 3 , 2 , a , r));
                            vaxToScreenTreatVaxHpv(d,v,s,aS,r) = toInd(allcomb(d , v , 2 , s , 1 , 4 , 2 , a , r));
                            noVaxToScreenTreatNonVaxHpv(d,v,h,aS,r) = toInd(allcomb(d , v , h , 2 , 1 , 3 , 2 , a , r));
                            vaxToScreenTreatNonVaxHpv(d,v,h,aS,r) = toInd(allcomb(d , v , h , 2 , 1 , 4 , 2 , a , r));
                            noVaxToScreenHyst(d,v,aS,r) = toInd(allcomb(d , v , 6 , 6 , 10 , 3 , 2 , a , r));
                            vaxToScreenHyst(d,v,aS,r) = toInd(allcomb(d , v , 6 , 6 , 10 , 4 , 2 , a , r));

                            % for people with cancer
                            noVaxToScreenCancerNoTreat(d,v,h,s,x,aS,r) = sort(toInd(allcomb(d , v , h , s , x+3 , 3 , 2 , a , r)));
                            noVaxToScreenCancerTreat(d,v,h,s,x,aS,r) = sort(toInd(allcomb(d , v , h , s , x+6 , 3 , 2 , a , r)));
                            noVaxToScreenCancerNegScreen(d,v,h,s,x,aS,r) = sort(toInd(allcomb(d , v , h , s , x , 3 , 2 , a , r)));
                            vaxToScreenCancerNoTreat(d,v,h,s,x,aS,r) = sort(toInd(allcomb(d , v , h , s , x+3 , 4 , 2 , a , r)));
                            vaxToScreenCancerTreat(d,v,h,s,x,aS,r) = sort(toInd(allcomb(d , v , h , s , x+6 , 4 , 2 , a , r)));
                            vaxToScreenCancerNegScreen(d,v,h,s,x,aS,r) = sort(toInd(allcomb(d , v , h , s , x , 4 , 2 , a , r)));
                        end
                    end    
                end
            end
        end

    end

    % Create indices for removing screening status as people age out of screened age groups
    noVaxScreen(:,aS) = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 3 , ... 
        2 , a+1 , 1 : risk));
    noVaxXscreen(:,aS) = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 , ... 
        2 , a+1 , 1 : risk));
    vaxScreen(:,aS) = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 4 , ... 
        2 , a+1 , 1 : risk));
    vaxXscreen(:,aS) = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 2 , ... 
        2 , a+1 , 1 : risk));
end

% Create indices for symptomatic CC detection
% Undetected local, regional, and distant population
for a = 1 : age
    for d = 1 : disease
        for v = 1 : viral
            for h = 1 : hpvVaxStates
                for s = 1 : hpvNonVaxStates
                    for x = 1 : 3 %only first 3 groups can be symptomatically detected 
                        for i = 1 : intervens
                            for r = 1 : risk
                                udPop(d,v,h,s,x,i,a,r) = sort(toInd(allcomb(d , v , h , s , x , i , 2 , a , r)));
                                udPopNoTreat(d,v,h,s,x,i,a,r) = sort(toInd(allcomb(d , v , h , s , x+3 , i , 2 , a , r)));
                                udPopTreat(d,v,h,s,x,i,a,r) = sort(toInd(allcomb(d , v , h , s , x+6 , i , 2 , a , r)));
                                udPopHyst(d,v,h,s,x,i,a,r) = sort(toInd(allcomb(d , v , 6 , 6 , 10 , i , 2 , a , r)));
                            end 
                        end 
                    end 
                end 
            end 
        end 
    end 
end 

%% Vaccination
lambdaMultVaxMat = zeros(age , 2);   % age-based vector for modifying lambda based on vaccination status

% No waning
lambdaMultVaxMat(min(vaxAge) : age , 1) = vaxEff;

% Waning
effPeriod = 20; % number of years that initial efficacy level is retained
wanePeriod = 20; % number of years over which initial efficacy level wanes
if waning 
    % Following a period (in years) where original efficacy is retained, 
    % specified by 'effPeriod' , linearly scale down vaccine efficacy 
    % to 0% over time period specificed by 'wanePeriod'
    % To make waning rate equal in all scenarios, the linear rate of 
    % waning is based on the least effective initial vaccine efficacy scenario.        
    kWane = vaxEff / round(wanePeriod / 5);     
    vaxInit = vaxEff;
    lambdaMultVaxMat(round(effPeriod / 5) + min(vaxAge) - 1 : age) = ...
        max(0 , linspace(vaxInit , ...
        vaxInit - kWane * (1 + age - (round(wanePeriod / 5) + min(vaxAge))) ,...
        age - (round(wanePeriod / 5) + min(vaxAge)) + 2)'); % ensures vaccine efficacy is >= 0
end
lambdaMultVax = 1 - lambdaMultVaxMat;

%% Simulation
% disp('Start up')
% profile on
% disp(' ')

% If starting from beginning
if ~ isfile([pwd , 'HHCoM_Results/' , pathModifier , '.mat'])
    
    % Initial Population 
    MpopStruc = riskDist(: , : , 1);
    FpopStruc = riskDist(: , : , 2);
    mPop = zeros(age , risk); % distribute initial population size by gender, age risk
    fPop = mPop;
    for i = 1 : age
        mPop(i , :) = MpopStruc(i, :).* mInit(i); % ./ (12*1.12); % scale population size back to an approximated level at the start year
        fPop(i , :) = FpopStruc(i, :).* fInit(i); % ./ (12*1.12);
    end
    initPop = zeros(dim);
    initPop(1 , 1 , 1 , 1 , 1 , 1 , 1 , : , :) = mPop; % HIV-, HPV Susceptible, no precancer, male
    initPop(1 , 1 , 1 , 1 , 1 , 1 , 2 , : , :) = fPop; % HIV-, HPV Susceptible, no precancer, female
    initPop_0 = initPop;
    % Assumes HPV starts before HIV in ages 15-44
    infected = initPop_0(1 , 1 , 1 , 1 , 1 , 1 , : , ...
        (15/max(1,fivYrAgeGrpsOn*5)+1) : (45/max(1,fivYrAgeGrpsOn*5)) , :) .* (0.2 * 0.9975); % initial HPV prevalence among ages 15-44 (sexually active) (HIV-)
    initPop(1 , 1 , 1 , 1 , 1 , 1 , : , (15/max(1,fivYrAgeGrpsOn*5)+1) : (45/max(1,fivYrAgeGrpsOn*5)) , :) = ...
        initPop_0(1 , 1 , 1 , 1 , 1 , 1 , 1 : gender , (15/max(1,fivYrAgeGrpsOn*5)+1) : (45/max(1,fivYrAgeGrpsOn*5)) , :) - infected; % moved from HPV susceptible
    %initPop(1 , 1 , 2 , 1 , 1 , 1 , : , (15/max(1,fivYrAgeGrpsOn*5)+1) : (45/max(1,fivYrAgeGrpsOn*5)) , :) = infected;
    initPop(1 , 1 , 2 , 1 , 1 , 1 , : , (15/max(1,fivYrAgeGrpsOn*5)+1) : (45/max(1,fivYrAgeGrpsOn*5)) , :) = 0.50 .* infected; % half moved to vaccine-type HPV+
    initPop(1 , 1 , 1 , 2 , 1 , 1 , : , (15/max(1,fivYrAgeGrpsOn*5)+1) : (45/max(1,fivYrAgeGrpsOn*5)) , :) = 0.45 .* infected; % moved to non-vaccine-type HPV+
    initPop(1 , 1 , 2 , 2 , 1 , 1 , : , (15/max(1,fivYrAgeGrpsOn*5)+1) : (45/max(1,fivYrAgeGrpsOn*5)) , :) = 0.05 .* infected; % moved to vaccine-type and non-vaccine-type HPV+  
    assert(~any(initPop(:) < 0) , 'Some compartments negative after seeding HPV infections.')
    popIn = reshape(initPop , prod(dim) , 1); % initial population to "seed" model
    
    % Directory to save results
    if ~ exist([pwd , '/HHCoM_Results/'])
        mkdir HHCoM_Results
    end
    
    % Initialize performance tracking vector
    % runtimes = zeros(size(s , 2) - 2 , 1);
    
    % Initialize time vector
    s = 1 : timeStep : years + 1 + timeStep;
    tVec = linspace(startYear , endYear , length(s) - 1);
    iStart = 2;
    disp(['Starting from beginning, year: ' , num2str(startYear + s(iStart) - 1) , ...
        ' and simulating to ' , num2str(endYear) , ' with ' , num2str(stepsPerYear) , ' steps per year.'])
    
    % Initialize result vectors
    popVec = spalloc(length(s) - 1 , prod(dim) , 10 ^ 8);
    popVec(1 , :) = popIn;
    deaths = zeros(length(s) - 1 , 1); %popVec; 
    newHiv = zeros(length(s) - 1 , hpvVaxStates , hpvNonVaxStates , endpoints , gender , age , risk);
    hivDeaths = zeros(length(s) - 1 , disease, gender , age);
    newHpvVax = zeros(length(s) - 1 , gender , disease , age , risk , intervens);
    newImmHpvVax = newHpvVax;
    newHpvNonVax = newHpvVax;
    newImmHpvNonVax = newHpvVax;
    newCC = zeros(length(s) - 1 , disease , age , hpvTypeGroups); % track by HPV type causal to CC
    % newCin1 = newCC;
    % newCin2 = newCC;
    % newCin3 = newCC;
    % ccDeath = newCC;
    ccDeath_treat = newCC; 
    ccDeath_untreat = newCC; 
    % ccDeath_treat_stage = zeros(length(s) - 1 , age , hpvTypeGroups , 6); 
    newScreen = zeros(length(s) - 1 , disease , viral , hpvVaxStates , hpvNonVaxStates , 3 , numScreenAge , 2);
%     newTreatImm = newScreen;
%     newTreatHpv = newScreen;
%     newTreatHyst = newScreen;
    menCirc = zeros(length(s) - 1 , 1);
    vaxdSchool = zeros(length(s) - 1 , 1);
    ccSymp = zeros(length(s) - 1 , 3 , age , 3); 
    ccTreat = ccSymp; 
    
    % ART
    import java.util.LinkedList
    artDistList = LinkedList();
    artDist = zeros(disease , viral , gender , age , risk); % initial distribution of inidividuals on ART = 0
    %artTreatTracker = zeros(length(s) - 1 , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , gender , age , risk);

% If continuing from checkpoint
elseif isfile([pwd , 'HHCoM_Results/' , pathModifier , '.mat'])
    % Initial Population 
    chckPntIn = load([pwd , '/HHCoM_Results/' , pathModifier]); % name for historical run input file 
    
    % Initialize time vector
    s = 1 : timeStep : years + 1 + timeStep;
    tVec = linspace(startYear , endYear , length(s) - 1);
    iStart = chckPntIn.i+1;
    disp(['Continuing from checkpoint, year: ' , num2str(startYear + s(iStart) - 1) , ...
        ' and simulating to ' , num2str(endYear) , ' with ' , num2str(stepsPerYear) , ' steps per year.'])
    
    % Initialize result vectors
    popVec = chckPntIn.popVec;
    deaths = chckPntIn.deaths; 
    newHiv = chckPntIn.newHiv;
    hivDeaths = chckPntIn.hivDeaths;
    newHpvVax = chckPntIn.newHpvVax;
    newImmHpvVax = chckPntIn.newImmHpvVax;
    newHpvNonVax = chckPntIn.newHpvNonVax;
    newImmHpvNonVax = chckPntIn.newImmHpvNonVax;
    newCC = chckPntIn.newCC; % track by HPV type causal to CC
    % newCin1 = chckPntIn.newCin1;
    % newCin2 = chckPntIn.newCin2;
    % newCin3 = chckPntIn.newCin3;
    % ccDeath = chckPntIn.ccDeath;
    ccDeath_treat = chckPntIn.ccDeath_treat; 
    ccDeath_untreat = chckPntIn.ccDeath_untreat; 
    % ccDeath_treat_stage = chckPntIn.ccDeath_treat_stage; 
    newScreen = chckPntIn.newScreen;
%     newTreatImm = chckPntIn.newTreatImm;
%     newTreatHpv = chckPntIn.newTreatHpv;
%     newTreatHyst = chckPntIn.newTreatHyst;
    menCirc = chckPntIn.menCirc;
    vaxdSchool = chckPntIn.vaxdSchool;
    ccSymp = chckPntIn.ccSymp; 
    ccTreat = chckPntIn.ccTreat;
    
    % ART
    import java.util.LinkedList
    artDistList = chckPntIn.artDistList;
    artDist = chckPntIn.artDist;
    %artTreatTracker = zeros(length(s) - 1 , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , gender , age , risk);
end

%% Main body of simulation
% disp(['Simulating period from ' num2str(startYear) ' to ' num2str(endYear) ...
%     ' with ' num2str(stepsPerYear), ' steps per year.'])
disp(' ')
disp('Simulation running...')
disp(' ')
% progressbar('Simulation Progress')

for i = iStart : length(s) - 1
    year = startYear + s(i) - 1
    tspan = [s(i) , s(i + 1)]; % evaluate diff eqs over one time interval
    popIn = popVec(i - 1 , :);
%   currStep = round(s(i) * stepsPerYear);
%   disp(['current step = ' num2str(startYear + s(i) - 1) ' ('...
%       num2str(length(s) - i) ' time steps remaining until year ' ...
%       num2str(endYear) ')'])
    
    % Add HIV index cases at hivStartYear
    if (hivOn && (year == hivStartYear))
        % Initialize HIV cases in population ages 15-29
        popIn_init = popIn;
        
        % Create indices
        % % fromNonHivNonHpv = sort(toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 1:gender , (15/max(1,fivYrAgeGrpsOn*5)+1) : (30/max(1,fivYrAgeGrpsOn*5)) , 1:risk))); 
        % % toHivNonHpv = sort(toInd(allcomb(4 , 2 , 1 , 1 , 1 , 1 , 1:gender , (15/max(1,fivYrAgeGrpsOn*5)+1) : (30/max(1,fivYrAgeGrpsOn*5)) , 1:risk)));
        % % fromNonHivHpv = sort(toInd(allcomb(1 , 1 , 2 : hpvVaxStates , 2 : hpvNonVaxStates , 1 : 3 , 1 , 1:gender , (15/max(1,fivYrAgeGrpsOn*5)+1) : (30/max(1,fivYrAgeGrpsOn*5)) , 1:risk))); 
        % % toHivHpv = sort(toInd(allcomb(4 , 2 , 2 : hpvVaxStates , 2 : hpvNonVaxStates , 1 : 3 , 1 , 1:gender , (15/max(1,fivYrAgeGrpsOn*5)+1) : (30/max(1,fivYrAgeGrpsOn*5)) , 1:risk)));
        fromNonHivAll1 = sort(toInd(allcomb(1 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 1 : gender , 1 : 3 , 1:risk)));
        fromNonHivAll2 = sort(toInd(allcomb(1 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 1 : gender , 4 : 7 , 1:risk)));
        fromNonHivAll3 = sort(toInd(allcomb(1 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 1 : gender , 8: 11, 1:risk)));
        fromNonHivAll4 = sort(toInd(allcomb(1 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 1 : gender , 12: age, 1:risk)));
        toHivAll1 = sort(toInd(allcomb(4 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 1 : gender , 1 : 3 , 1:risk)));
        toHivAll2 = sort(toInd(allcomb(4 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 1 : gender , 4 : 7 , 1:risk)));
        toHivAll3 = sort(toInd(allcomb(4 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 1 : gender , 8 : 11 , 1:risk)));
        toHivAll4 = sort(toInd(allcomb(4 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 1 : gender , 12 : age , 1:risk)));

        % Distribute HIV infections 
        % % (HPV-)        
        % % popIn(fromNonHivNonHpv) = (1 - 0.002) .* popIn_init(fromNonHivNonHpv);    % reduce non-HIV infected
        % % popIn(toHivNonHpv) = (0.002) .* popIn_init(fromNonHivNonHpv);    % increase HIV infected ( male/female, age groups 4-6) (% prevalence)
        % % (HPV+)
        % % popIn(fromNonHivHpv) = (1 - 0.001) .* popIn_init(fromNonHivHpv);    % reduce non-HIV infected
        % % popIn(toHivHpv) = (0.001) .* popIn_init(fromNonHivHpv);    % increase HIV infected ( male/female, age groups 4-6) (% prevalence)
        popIn(fromNonHivAll1) = (1 - 0.0002) .* popIn_init(fromNonHivAll1);    % reduce non-HIV infected
        popIn(toHivAll1) = (0.0002) .* popIn_init(fromNonHivAll1);    % increase HIV infected ( male/female, age groups 4-6) (% prevalence)
        popIn(fromNonHivAll2) = (1 - 0.005) .* popIn_init(fromNonHivAll2);    % reduce non-HIV infected
        popIn(toHivAll2) = (0.005) .* popIn_init(fromNonHivAll2);    % increase HIV infected ( male/female, age groups 4-6) (% prevalence)
        popIn(fromNonHivAll3) = (1 - 0.006) .* popIn_init(fromNonHivAll3);    % reduce non-HIV infected
        popIn(toHivAll3) = (0.006) .* popIn_init(fromNonHivAll3);    % increase HIV infected ( male/female, age groups 4-6) (% prevalence)
        popIn(fromNonHivAll4) = (1 - 0.0002) .* popIn_init(fromNonHivAll4);    % reduce non-HIV infected
        popIn(toHivAll4) = (0.0002) .* popIn_init(fromNonHivAll4);
    end

    if hpvOn
        % HPV NATURAL HISTORY
        % Progression and clearance of HPV
        % Progression and regression of precancerous lesions
        % Development and progression of cervical cancer
        % Cervical cancer-associated mortality by stage and HIV status/CD4 count
        [~ , pop , newCC(i , : , : , :) , ccDeath_untreat(i , : , : , :), ccDeath_treat(i , : , : , :)] ...
            = ode4xtra(@(t , pop) ...
            hpvCCNH(t , pop , hpv_hivClear , rImmuneHiv , c3c2Mults , c2c1Mults , muCC , muCC_ud , muCC_d , ...
            normalhpvVaxInds , immunehpvVaxInds , infhpvVaxInds , normalhpvNonVaxInds , ...
            immunehpvNonVaxInds , infhpvNonVaxInds , cin3hpvVaxIndsFrom , ccLochpvVaxIndsTo , ...
            ccLochpvVaxIndsFrom , ccReghpvVaxInds , ccDisthpvVaxInds , ...
            cin3hpvNonVaxIndsFrom , ccLochpvNonVaxIndsTo , ccLochpvNonVaxIndsFrom , ...
            ccReghpvNonVaxInds , ccDisthpvNonVaxInds , cin1hpvVaxInds , ...
            cin2hpvVaxInds , cin3hpvVaxInds , cin1hpvNonVaxInds , ...
            cin2hpvNonVaxInds , cin3hpvNonVaxInds , kInf_Cin1 , kCin1_Cin2 , kCin2_Cin3 , ...
            kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , kCin1_Inf , rNormal_Inf , ...
            rImmune , fImm , kRL , kDR , maleHpvClearMult , disease , age , hpvVaxStates , hpvNonVaxStates , hpvTypeGroups , ccLochpvVaxIndsFrom_treat , ...
            ccReghpvVaxInds_treat , ccDisthpvVaxInds_treat) , tspan , popIn);
        popIn = pop(end , :); % for next module
        if any(pop(end , :) <  0)
            disp('After hpv')
            break
        end
        pop = pop(end , :); % next module reads in pop, not popInd
        
        if (year >= hpvScreenStartYear)
            % CERVICAL CANCER SCREENING AND TREATMENT
            % Screening
            % Treatment
            [dPop , newScreen(i , : , : , : , :, : , : , :), ccTreat(i, : , : , :)]   ...
                = hpvScreen(pop , ...
                    disease , viral , age , hpvVaxStates , hpvNonVaxStates , intervens , endpoints , risk , ...
                    screenYrs , screenAlgs , year , stepsPerYear , screenAgeAll , screenAgeS , ...
                    noVaxNoScreen , noVaxToScreen , vaxNoScreen , vaxToScreen , noVaxToScreenTreatImm , ...
                    vaxToScreenTreatImm , noVaxToScreenTreatHpv , vaxToScreenTreatHpv , ...
                    noVaxToScreenTreatVaxHpv , vaxToScreenTreatVaxHpv , noVaxToScreenTreatNonVaxHpv , ...
                    vaxToScreenTreatNonVaxHpv , noVaxToScreenHyst , vaxToScreenHyst , numScreenAge , ageMultsComb , ...
                    noVaxToScreenCancerNoTreat , noVaxToScreenCancerTreat , ...
                    vaxToScreenCancerNoTreat , vaxToScreenCancerTreat , hystMult , kSymp , udPop , udPopNoTreat , udPopTreat , udPopHyst , ...
                    vaxToScreenCancerNegScreen , noVaxToScreenCancerNegScreen);
            pop(end , :) = pop(end , :) + dPop;
            popIn = pop(end , :); % for next module
            if any(pop(end , :) <  0)
                disp('After hpv screen')
                break
            end
        end
        pop = pop(end, :); 

        % SYMPTOMATIC CC DETECTION IN A NON-SCREENING YEAR
%         if (year < hpvScreenStartYear)
            [dPop , ccSymp(i,:,:,:)] = symptomaticDetection(pop , ...
                year , hpvScreenStartYear , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , risk , intervens , age , ...
                screenAlgs , kSymp, hystMult, udPop, udPopNoTreat, udPopTreat, udPopHyst, ccReghpvVaxInds); 
            pop(end , :) = pop(end , :) + dPop(end , :); 
            popIn = pop(end , :); % for next module
            if any(pop(end , :) < 0)
                disp('After symptomatic detection')
                break
            end 
%         end 
        pop = pop(end, :); 
    end
    
    % HPV AND HIV TRANSMISSION
    % Heterosexual mixing by gender, age, and risk group
    % Partnership adjustment
    % HPV infection by type
    % HIV infection and protection provided by condoms, circumcision, and ART
    [~ , pop , newHpvVax(i , : , : , : , : , :) , newImmHpvVax(i , : , : , : , : , :) , ...
        newHpvNonVax(i , : , : , : , : , :) , newImmHpvNonVax , newHiv(i , : , : , : , : , : , :)] = ...
        ode4xtra(@(t , pop) mixInfect(t , pop , ...
        stepsPerYear , year , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , intervens , gender , ...
        age , risk , fivYrAgeGrpsOn , hpvTypeGroups , ageSexDebut , gar , epsA_vec , epsR_vec , yr , ...
        partnersM , partnersF , partnersMmult,...
        beta_hpvVax_mod , beta_hpvNonVax_mod , vaxInds , nonVInds , ...
        lambdaMultImm , lambdaMultVax , artHpvMult , hpv_hivMult , ...
        hpvVaxSus , hpvVaxImm , hpvVaxInf , hpvNonVaxSus , hpvNonVaxImm , hpvNonVaxInf , ...
        circProtect , condProtect , condUse , betaHIV_mod , hiv_hpvMult, ...
        d_partnersMmult,  ...
        hivSus , toHiv , hivCurr) , tspan , popIn);
    popIn = pop(end , :); % for next module
    if any(pop(end , :) < 0)
        disp('After mixInfect')
        break
    end
    
    % HIV NATURAL HISTORY
    % CD4 progression
    % Viral load progression
    % ART initiation, dicontinuation, and scale-up by CD4 count
    % HIV-associated mortality
    if (hivOn && (year >= hivStartYear))
        [~ , pop , hivDeaths(i , :, : , :) , artTreat] =...
            ode4xtra(@(t , pop) hivNH(t , pop , vlAdvancer , muHIV , dMue , mue3 , mue4 , artDist , ...
            kCD4 ,  artYr_vec , artM_vec , artF_vec , minLim , maxLim , disease , viral , ...
            hpvVaxStates , hpvNonVaxStates , endpoints , gender , age , risk , ...
            ageSexDebut , hivInds , stepsPerYear , year) , tspan , popIn);
        popIn = pop(end , :);
        %artTreatTracker(i , : , : , : , : , : , : , :  ,:) = artTreat;
        artDistList.add(sum(sum(sum(artTreat , 3) , 4) , 5));
        if artDistList.size() >= stepsPerYear * 2
            artDistList.remove(); % remove CD4 and VL distribution info for people initiating ART more than 2 years ago
        end
        artDist = calcDist(artDistList , disease , viral , gender , age , ...
            risk);
        if any(pop(end , :) < 0)
            disp('After hiv')
            break
        end
    end
    
    % DEMOGRAPHY
    % Births
    % Mother-to-child HIV transmission
    % Aging and risk-group redistribution
    % Natural deaths
    [~ , pop , deaths(i , :)] = ode4xtra(@(t , pop) ...
        bornAgeDieRisk(t , pop , year , ...
        gender , age , fivYrAgeGrpsOn , fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
        hivFertPosBirth2 , hivFertNegBirth2 , fertMat3 , hivFertPosBirth3 , hivFertNegBirth3 , ...
        fertMat4 , hivFertPosBirth4 , hivFertNegBirth4 , ...
        dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ...
        dFertPos3 , dFertNeg3  , dFertMat3,  ...
        deathMat , deathMat2 , deathMat3 , deathMat4 , deathMat5, ...
        dDeathMat , dDeathMat2 , dDeathMat3, dDeathMat4,...
        MTCTRate , ageInd , riskAdj, d_riskAdj, riskInd , riskDist ,......
        stepsPerYear , currYear , agesComb , noVaxScreen , noVaxXscreen , ...
        vaxScreen , vaxXscreen , hpvScreenStartYear) , tspan , popIn);
    popIn = pop(end , :);
    if any(pop(end , :) < 0)
        disp('After bornAgeDieRisk')
        break
    end 
    
    % VOLUNTARY MALE MEDICAL CIRCUMCISION
    % Scale-up of VMMC by age
    if (year >= circStartYear)
        [dPop , menCirc(i , :)] = vmmc(popIn , circStartYear , circNatStartYear , ...
            vmmcYr_vec , vmmc_vec , circ_aVec , hivNegNonVMMCinds , hivNegVMMCinds , ...
            ageSexDebut , year);
        pop(end , :) = pop(end , :) + dPop;
        popIn = pop(end , :);
        if any(pop(end , :) < 0)
            disp('After vmmc')
            break
        end
    end
    
    if ((year >= vaxStartYear)) %&& (vaxRate > 0))
        % HPV VACCINATION
        % School-based vaccination regimen
        [dPop , vaxdSchool(i , :)] = hpvVaxSchool(popIn , disease , viral , risk , ...
            hpvVaxStates , hpvNonVaxStates , endpoints , intervens , vaxG , vaxAge , ...
            vaxRate_vec , toInd , vaxYrs , year , stepsPerYear , gradScaleUp); 
        pop(end , :) = pop(end , :) + dPop;
        if any(pop(end , :) < 0)
            disp('After hpvVaxSchool')
            break
        end
    end

    % add results to population vector
    popVec(i , :) = pop(end , :)';
    % runtimes(i) = toc;
    % progressbar(i/(length(s) - 1))
    
%     if rem(year , 25) == 0.0
%         savdir = [pwd , '/HHCoM_Results/'];
%         save(fullfile(savdir , pathModifier, '') , 'fivYrAgeGrpsOn' , 'tVec' ,  'popVec' , 'newHiv' , ...
%             'newHpvVax' , 'newImmHpvVax' , 'newHpvNonVax' , 'newImmHpvNonVax' , ...
%             'hivDeaths' , 'deaths' ,  'menCirc' , 'vaxdSchool' , ...
%             'newScreen' , 'ccDeath_treat', 'ccDeath_untreat' ,  ... % 'newTreatImm' , 'newTreatHpv' , 'newTreatHyst' , ...
%             'newCC' , 'artDist' , 'artDistList' , ... % 'artTreatTracker' , ...
%             'ccSymp' , 'ccTreat' , ...
%             'startYear' , 'endYear' , 'i' , '-v7.3');
%     end

    % Checking for year 1999 to check stage distribution
    if year >= 1950 && year <= 1950.1
        local = sum(ccSymp(i, 1, 1:end, 1:end), "all"); 
        regional = sum(ccSymp(i, 2, 1:end, 1:end), "all"); 
        distant = sum(ccSymp(i, 3, 1:end, 1:end), "all"); 
        total = local+regional+distant; 
        stageDist = [local/total regional/total distant/total]; 
    end 

    % Checking for year 2022 to check the total number of cancers detected
    if year >= 1950 && year <= 1950.1
        local = sum(ccSymp(i, 1, 1:end, 1:end), "all") + sum(ccTreat(i, 1, 1:end, 1:end), "all"); 
        regional = sum(ccSymp(i, 2, 1:end, 1:end), "all") + sum(ccTreat(i, 2, 1:end, 1:end), "all"); 
        distant = sum(ccSymp(i, 3, 1:end, 1:end), "all") + sum(ccTreat(i, 3, 1:end, 1:end), "all"); 
        total = local+regional+distant; 
        numDxCC = [local regional distant total];
    end 
    
    disp(['Reached year ' num2str(year)])

end
popLast = sparse(popVec(end-1 , :));
disp(['Reached year ' num2str(endYear)])
popVec = sparse(popVec); % compress population vectors

%% Save results
savdir = [pwd , '/HHCoM_Results/'];
save(fullfile(savdir , pathModifier, '') , 'stageDist' , 'numDxCC', 'kSymp', 'kRL', 'kDR');

disp(' ')
disp('Simulation complete.')
toc
% profile viewer

%% Calculate summed log-likelihood
% if calibBool    
%     negSumLogL = likeFun(popVec , newCC , cinPos2007_dObs , cin1_2010_dObs ,...
%         cin2_2010_dObs, hpv_hiv_dObs , hpv_hivNeg_dObs , hivPrevM_dObs , hivPrevF_dObs , ...
%         hivPrevAll_dObs, hpv_all_dObs , hpv_hiv2009_dObs , cc_dist_dObs , cin3_dist_dObs , ...
%         cin1_dist_dObs , hpv_dist_dObs , ccInc2012_dObs, popAgeDist_dObs , totPopSize_dObs , ...
%         toInd , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , intervens , ...
%         age , gender , risk , startYear , stepsPerYear , annlz)
% 
%     if negSumLogL < -190000.00
%         delete([savdir , pathModifier , '.mat']);
%     end
%     
%     if isnan(negSumLogL)
%         negSumLogL = -10000000.00;
%     end
%     
% else
%     negSumLogL = likeFun(popVec , newCC , cinPos2007_dObs , cin1_2010_dObs ,...
%         cin2_2010_dObs, hpv_hiv_dObs , hpv_hivNeg_dObs , hivPrevM_dObs , hivPrevF_dObs , ...
%         hivPrevAll_dObs, hpv_all_dObs , hpv_hiv2009_dObs , cc_dist_dObs , cin3_dist_dObs , ...
%         cin1_dist_dObs , hpv_dist_dObs , ccInc2012_dObs, popAgeDist_dObs , totPopSize_dObs , ...
%         toInd , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , intervens , ...
%         age , gender , risk , startYear , stepsPerYear , annlz)
% end

%% Runtimes
% figure()
% plot(1 : size(runtimes , 1) , runtimes)
% xlabel('Step'); ylabel('Time(s)')
% title('Runtimes')
% avgRuntime = mean(runtimes); % seconds
% stdRuntime = std(runtimes); % seconds
% disp(['Total runtime: ' , num2str(sum(runtimes) / 3600) , ' hrs' , ' (' , num2str(sum(runtimes) / 60) , ' mins)']);
% disp(['Average runtime per step: ' , num2str(avgRuntime / 60) , ' mins (' , num2str(avgRuntime) , ' secs)']);
% disp(['Standard deviation: ' , num2str(stdRuntime / 60) , ' mins (' , num2str(stdRuntime) , ' secs)']);
% figure()
% h = histogram(runtimes);
% title('Runtimes')
% ylabel('Frequency')
% xlabel('Times (s)')

%% Show results
% negSumLogL
% showResults(pathModifier)
