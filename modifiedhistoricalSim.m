% Historical module
% Runs simulation over the time period and time step specified by the user.
% Modified for the symptomatic rate calibration

function [negSumLogL] = modifiedhistoricalSim(calibBool , pIdx , paramsSub , paramSet , paramSetIdx , tstep_abc , date , sympParams_in , sympRun)    % input variables when using a calibration parameter set
% historicalSim(0 , [] , [] , [] , [] , 0 , '19May20')    % input variables when running from command window using hand-calibrated, hard-coded parameter values
% Note: if you hard-code the "pathModifier" file output name variable below, then the date, paramSetIdx, and tstep_abc input values here are just dummy values and unused

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
% pathModifier = ['toNow_' , date , '_2v57BaseVax_spCytoScreen_shortName_noVMMChpv_discontFxd_screenCovFxd_hivInt2017_' , num2str(tstep_abc) , '_' , num2str(paramSetIdx)]; % ***SET ME***: name for historical run output file 
pathModifier = ['toNow_', date, '_testingTreatment_changeSympProb_1925to1945_', num2str(sympRun), '_' num2str(paramSetIdx)]; 
%pathModifier = ['toNow_' , date , '_baseVax057_baseScreen_baseVMMC_DoART_S3_' , num2str(tstep_abc) , '_' , num2str(paramSetIdx)]; % ***SET ME***: name for historical run output file 
%pathModifier = 'toNow_21Feb20_testMuART_1925Start_decBkrndMort';

% AGE GROUPS
fivYrAgeGrpsOn = 1; % choose whether to use 5-year (fivYrAgeGrpsOn=1) or 1-year age groups (fivYrAgeGrpsOn=0)

% SCREENING
% Instructions: Choose either the original baseline screening algorithm or the screening paper cytology algorithm
%   Screening coverage, HIV groups, and ages are set automatically below 
screenAlgorithm = 3; % ***SET ME***: screening algorithm to use (1 for baseline, 3 for spCyto)

% VACCINATION
% Instructions: If you want no historical vaccination/ no vaccination in your baseline scenario, set baseline vaccine coverage to zero. 
%   Otherwise, set it to 0.86 (the historical coverage level for 9 year-olds) times an adjustment for the bivalent vaccine. We are 
%   adjusting coverage rather than efficacy because we don't track vaccination by vaccine type and age cohorts can have mixed vaccine
%   types, particularly if catch-up vaccination is applied in future years.

% Common parameters
vaxEff = 1.0; % actually bivalent vaccine, but to avoid adding additional compartments, we use nonavalent vaccine and then reduce coverage
rVaxWane = 0.0; % rate of waning vaccine immunity

%Parameters for school-based vaccination regimen  % ***SET ME***: coverage for baseline vaccination of 9-year-old girls
vaxAge = 2;
vaxRate = 0.57*(0.7/0.9); %0.86*(0.7/0.9);    % (9 year-old coverage * bivalent vaccine efficacy adjustment (0.7/0.9 proportion of cancers prevented)); last dose, first dose pilot
vaxG = 2;   % indices of genders to vaccinate (1 or 2 or 1,2); set stepsPerYear=8 in loadUp2.m if including vaccination of boys 

% ART + VIRAL SUPPRESSION & VMMC COVERAGE
% Instructions: In loadUp2.m, go to the section titled "Save intervention
%   parameters." In the sub-sections labeled "ART+VS coverage" and "VMMC
%   coverage", select your desired scale-up assumptions    % ***SET ME***: ART & VMMC scale-up assumptions

%% Save pre-loaded parameters and pre-calculated indices and matrices
[stepsPerYear , timeStep , startYear , currYear , endYear , ...
    years , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , ...
    intervens , gender , age , risk , hpvTypeGroups , dim , k , toInd , ...
    annlz , ...
    ageSexDebut , mInit , fInit , partnersM , partnersF , maleActs , ...
    femaleActs , riskDist , fertility , fertility2 , fertility3 , fertility4 , ...
    mue , mue2 , mue3 , mue4 , epsA_vec , epsR_vec , ...
    yr , ...
    hivOn , betaHIV_mod , muHIV , kCD4 , ...
    hpvOn , beta_hpvVax_mod , beta_hpvNonVax_mod , fImm , rImmune , ...
    kCin1_Inf , kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , rNormal_Inf , kInf_Cin1 , ...
    kCin1_Cin2 , kCin2_Cin3 , lambdaMultImm , hpv_hivClear , rImmuneHiv , ...
    c3c2Mults , c2c1Mults , c2c3Mults , c1c2Mults , muCC , muCC_ud, muCC_d, kRL , kDR , artHpvMult , ...
    hpv_hivMult , maleHpvClearMult , ...
    condUse , screenYrs , hpvScreenStartYear , ...
    artYr , maxRateM , maxRateF , ...
    artYr_vec , artM_vec , artF_vec , minLim , maxLim , ...
    circ_aVec , vmmcYr_vec , vmmc_vec , vmmcYr , vmmcRate , ...
    hivStartYear , circStartYear , circNatStartYear , vaxStartYear , ...
    baseline , who , spCyto , spHpvDna , spGentyp , spAve , spHpvAve , ...
    circProtect , condProtect , MTCTRate , hyst , ...
    OMEGA , ...
    ccInc2012_dObs , ccInc2018_dObs , cc_dist_dObs , cin3_dist_dObs , ...
    cin1_dist_dObs , hpv_dist_dObs , cinPos2002_dObs , cinNeg2002_dObs , ...
    cinPos2015_dObs , cinNeg2015_dObs , hpv_hiv_dObs , hpv_hivNeg_dObs , ...
    hpv_hivM2008_dObs , hpv_hivMNeg2008_dObs , hivPrevM_dObs , hivPrevF_dObs , ...
    popAgeDist_dObs , totPopSize_dObs , ...
    hivCurr , ...
    gar , hivSus , hpvVaxSus , hpvVaxImm , hpvNonVaxSus , hpvNonVaxImm , ...
    toHiv , vaxInds , nonVInds , hpvVaxInf , hpvNonVaxInf , ...
    hivInds , ...
    cin3hpvVaxIndsFrom , ccLochpvVaxIndsTo , ccLochpvVaxIndsFrom , ...
    ccReghpvVaxInds , ccDisthpvVaxInds , cin3hpvNonVaxIndsFrom , ...
    ccLochpvNonVaxIndsTo , ccLochpvNonVaxIndsFrom , ccReghpvNonVaxInds , ...
    ccDisthpvNonVaxInds , cin1hpvVaxInds , cin2hpvVaxInds , cin3hpvVaxInds , ...
    cin1hpvNonVaxInds , cin2hpvNonVaxInds , cin3hpvNonVaxInds , normalhpvVaxInds , ...
    immunehpvVaxInds , infhpvVaxInds , normalhpvNonVaxInds , immunehpvNonVaxInds , ...
    infhpvNonVaxInds , fromVaxNoScrnInds , fromVaxScrnInds , toNonVaxNoScrnInds , ...
    toNonVaxScrnInds , ageInd , riskInd , ...
    kSymp , hystMult , ...
    hivNegNonVMMCinds , hivNegVMMCinds , ...
    vlAdvancer , ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , fertMat3 , hivFertPosBirth3 , hivFertNegBirth3 , ...
    fertMat4 , hivFertPosBirth4 , hivFertNegBirth4 , ...
    dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ...
    dFertPos3 , dFertNeg3 , dFertMat3 , deathMat , deathMat2 , deathMat3 , deathMat4 , ...
    dDeathMat , dDeathMat2 , dDeathMat3 , dMue, ...
    ccLochpvVaxIndsFrom_treat , ...
    ccReghpvVaxInds_treat , ccDisthpvVaxInds_treat] = loadUp2(fivYrAgeGrpsOn , calibBool , pIdx , paramsSub , paramSet);

%% Modifying the rate of symptomatic detection

kSymp = sympParams_in; 

%% Screening
if (screenAlgorithm == 1)
    % Baseline screening algorithm
    screenAlgs = baseline;
elseif (screenAlgorithm == 3)
    % Screening paper cytology algorithm
    screenAlgs = spCyto;
end
screenAlgs.genTypBool = 0;
screenAlgs.screenHivGrps = {[1:disease]};
screenAlgs.screenAge = {8};
screenAlgs.screenCover = [0.0; 0.18; 0.48; 0.48; 0.48; 0.48; 0.48]; % Coverage over time (Years: [2000; 2003; 2016; currYear; 2023; 2030; 2045])
screenAlgs.screenCover_vec = cell(size(screenYrs , 1) - 1, 1); % save data over time interval in a cell array
for i = 1 : size(screenYrs , 1) - 1          % interpolate dnaTestCover values at steps within period
    period = [screenYrs(i) , screenYrs(i + 1)];
    screenAlgs.screenCover_vec{i} = interp1(period , screenAlgs.screenCover(i : i + 1 , 1) , ...
        screenYrs(i) : timeStep : screenYrs(i + 1));
end

% Create screening indices
numScreenAge = length(screenAlgs.screenAge);
agesComb = screenAlgs.screenAge{1};
screenAgeAll = zeros(disease , viral , numScreenAge , risk , hpvVaxStates*hpvNonVaxStates*3*intervens);
screenAgeS = zeros(disease , viral , numScreenAge , risk , hpvVaxStates*hpvNonVaxStates*3*2);
noVaxNoScreen = zeros(disease , viral , hpvVaxStates , hpvNonVaxStates , 3 , numScreenAge , risk);
noVaxToScreen = noVaxNoScreen;
vaxNoScreen = noVaxNoScreen;
vaxToScreen = noVaxNoScreen;
noVaxToScreenTreatImm = zeros(disease , viral , numScreenAge , risk);
vaxToScreenTreatImm = noVaxToScreenTreatImm;
noVaxToScreenTreatHpv = noVaxToScreenTreatImm;
vaxToScreenTreatHpv = noVaxToScreenTreatImm;
noVaxToScreenHyst = noVaxToScreenTreatImm;
vaxToScreenHyst = noVaxToScreenTreatImm;
noVaxToScreenTreatImmVaxHpv = zeros(disease , viral , hpvNonVaxStates , numScreenAge , risk);
vaxToScreenTreatImmVaxHpv = noVaxToScreenTreatImmVaxHpv;
noVaxToScreenTreatImmNonVaxHpv = zeros(disease , viral , hpvVaxStates , numScreenAge , risk);
vaxToScreenTreatImmNonVaxHpv = noVaxToScreenTreatImmNonVaxHpv;
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
vaxToScreenCancerNoTreat = noVaxToScreenCancerNoTreat;
vaxToScreenCancerTreat = noVaxToScreenCancerNoTreat;
noVaxToScreenCancerNegScreen = noVaxToScreenCancerNoTreat; 
vaxToScreenCancerNegScreen = noVaxToScreenCancerNoTreat; 
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
                    for xS = 1 : 3 %only first 3 groups can get screened
                        x = xS; 
                        for r = 1 : risk
                            screenAgeAll(d,v,aS,r,:) = toInd(allcomb(d , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : 3 , 1 : intervens , 2 , a , r)); 
                            screenAgeS(d,v,aS,r,:) = toInd(allcomb(d , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : 3 , 3 : intervens , 2 , a , r));

                            noVaxNoScreen(d,v,h,s,xS,aS,r) = sort(toInd(allcomb(d , v , h , s , x , 1 , 2 , a , r)));
                            noVaxToScreen(d,v,h,s,xS,aS,r) = sort(toInd(allcomb(d , v , h , s , x , 3 , 2 , a , r)));
                            vaxNoScreen(d,v,h,s,xS,aS,r) = sort(toInd(allcomb(d , v , h , s , x , 2 , 2 , a , r)));
                            vaxToScreen(d,v,h,s,xS,aS,r) = sort(toInd(allcomb(d , v , h , s , x , 4 , 2 , a , r)));

                            % CC TREATMENT EDIT - keeping x=1 because that
                            % only undiagnosed will likely be treated
                            noVaxToScreenTreatImm(d,v,aS,r) = toInd(allcomb(d , v , 7 , 7 , 1 , 3 , 2 , a , r));
                            vaxToScreenTreatImm(d,v,r) = toInd(allcomb(d , v , 7 , 7 , 1 , 4 , 2 , a , r));
                            noVaxToScreenTreatImmVaxHpv(d,v,s,aS,r) = toInd(allcomb(d , v , 7 , s , 1 , 3 , 2 , a , r));
                            vaxToScreenTreatImmVaxHpv(d,v,s,aS,r) = toInd(allcomb(d , v , 7 , s , 1 , 4 , 2 , a , r));
                            noVaxToScreenTreatImmNonVaxHpv(d,v,h,aS,r) = toInd(allcomb(d , v , h , 7 , 1 , 3 , 2 , a , r));
                            vaxToScreenTreatImmNonVaxHpv(d,v,h,aS,r) = toInd(allcomb(d , v , h , 7 , 1 , 4 , 2 , a , r));
                            noVaxToScreenTreatHpv(d,v,aS,r) = toInd(allcomb(d , v , 2 , 2 , 1 , 3 , 2 , a , r));
                            vaxToScreenTreatHpv(d,v,aS,r) = toInd(allcomb(d , v , 2 , 2 , 1 , 4 , 2 , a , r));
                            noVaxToScreenTreatVaxHpv(d,v,s,aS,r) = toInd(allcomb(d , v , 2 , s , 1 , 3 , 2 , a , r));
                            vaxToScreenTreatVaxHpv(d,v,s,aS,r) = toInd(allcomb(d , v , 2 , s , 1 , 4 , 2 , a , r));
                            noVaxToScreenTreatNonVaxHpv(d,v,h,aS,r) = toInd(allcomb(d , v , h , 2 , 1 , 3 , 2 , a , r));
                            vaxToScreenTreatNonVaxHpv(d,v,h,r) = toInd(allcomb(d , v , h , 2 , 1 , 4 , 2 , a , r));
                            noVaxToScreenHyst(d,v,aS,r) = toInd(allcomb(d , v , 6 , 6 , 10 , 3 , 2 , a , r));
                            vaxToScreenHyst(d,v,aS,r) = toInd(allcomb(d , v , 6 , 6 , 10 , 4 , 2 , a , r));
                            
                            % for people with cancer
                            noVaxToScreenCancerNoTreat(d,v,h,s,xS,aS,r) = sort(toInd(allcomb(d , v , h , s , x+3 , 3 , 2 , a , r)));
                            noVaxToScreenCancerTreat(d,v,h,s,xS,aS,r) = sort(toInd(allcomb(d , v , h , s , x+6 , 3 , 2 , a , r)));
                            noVaxToScreenCancerNegScreen(d,v,h,s,xS,aS,r) = sort(toInd(allcomb(d , v , h , s , x , 3 , 2 , a , r)));
                            vaxToScreenCancerNoTreat(d,v,h,s,xS,aS,r) = sort(toInd(allcomb(d , v , h , s , x+3 , 4 , 2 , a , r)));
                            vaxToScreenCancerTreat(d,v,h,s,xS,aS,r) = sort(toInd(allcomb(d , v , h , s , x+6 , 4 , 2 , a , r)));
                            vaxToScreenCancerNegScreen(d,v,h,s,xS,aS,r) = sort(toInd(allcomb(d , v , h , s , x , 4 , 2 , a , r)));

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
lambdaMultVaxMat(min(vaxAge) : age , 1) = vaxEff;
lambdaMultVax = 1 - lambdaMultVaxMat;

%% Simulation
% disp('Start up')
% profile on
% disp(' ')

% If starting from beginning
if ~ isfile([pwd , '/HHCoM_Results/' , pathModifier , '.mat'])
    
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
    hivDeaths = zeros(length(s) - 1 , disease , gender , age);
    newHpvVax = zeros(length(s) - 1 , gender , disease , age , risk , intervens);
    newImmHpvVax = newHpvVax;
    newHpvNonVax = newHpvVax;
    newImmHpvNonVax = newHpvVax;
    newCC = zeros(length(s) - 1 , disease , age , hpvTypeGroups); % track by HPV type causal to CC
    % newCin1 = newCC;
    % newCin2 = newCC;
    % newCin3 = newCC;
    ccDeath = newCC;
    ccDeath_treat = newCC; 
    ccDeath_untreat = newCC; 
    ccDeath_treat_stage = zeros(length(s) - 1 , disease , age , hpvTypeGroups , 6); 
    newScreen = zeros(length(s) - 1 , disease , hpvVaxStates , hpvNonVaxStates , 3 , numScreenAge , 2);
%     newUndx = zeros(length(s) - 1 , disease , hpvVaxStates , hpvNonVaxStates , endpoints , age , 2); 
%     newDxUntreat = zeros(length(s) - 1 , disease , hpvVaxStates , hpvNonVaxStates , endpoints , age , 2); 
%     newTreat = zeros(length(s) - 1 , disease , hpvVaxStates , hpvNonVaxStates , endpoints , age , 2); % TODO: should i limit to just screening ages? 
    % newTreatImm = newScreen;
    % newTreatHpv = newScreen;
    % newTreatHyst = newScreen;
    % ccScreenNoTreat = newScreen; 
    % ccScreenTreat = newScreen; 
    % ccTreatHyst = newScreen; 
    menCirc = zeros(length(s) - 1 , 1);
    vaxdSchool = zeros(length(s) - 1 , 1);
    ccSymp = zeros(length(s) - 1 , disease , hpvVaxStates , hpvNonVaxStates , 3 , intervens , age , 3); 
    ccTreat = ccSymp; 
    toScreenMult_collect = zeros(length(s)-1, disease, viral, hpvVaxStates, hpvNonVaxStates, endpoints, age, risk);
    toScreenNoTreat_collect = toScreenMult_collect;
    toScreenNeg_collect = toScreenMult_collect; 
    toScreenTreat_collect = toScreenMult_collect; 
    toScreenTreatHystMult_collect = toScreenMult_collect; 
    
    % ART
    import java.util.LinkedList
    artDistList = LinkedList();
    artDist = zeros(disease , viral , gender , age , risk); % initial distribution of inidividuals on ART = 0
    artTreatTracker = zeros(length(s) - 1 , disease , viral , gender , age , risk); %zeros(length(s) - 1 , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , gender , age , risk);
    artDiscont = zeros(length(s) - 1 , disease , viral , gender , age , risk);
    
% If continuing from checkpoint
elseif isfile([pwd , '/HHCoM_Results/' , pathModifier , '.mat'])
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
    ccDeath = chckPntIn.ccDeath;
    ccDeath_treat = chckPntIn.ccDeath_treat; 
    ccDeath_untreat = chckPntIn.ccDeath_untreat; 
    ccDeath_treat_stage = chckPntIn.ccDeath_treat_stage; 
    newScreen = chckPntIn.newScreen;
%     newUndx = chckPntIn.newUndx; 
%     newUndxUntreat = chckPntIn.newDxUntreat; 
%     newTreat = chckPntIn.newTreat; 
    % newTreatImm = chckPntIn.newTreatImm;
    % newTreatHpv = chckPntIn.newTreatHpv;
    % newTreatHyst = chckPntIn.newTreatHyst;
    menCirc = chckPntIn.menCirc;
    vaxdSchool = chckPntIn.vaxdSchool;
    ccSymp = chckPntIn.ccSymp; 
    ccTreat = chckPntIn.ccTreat;
    toScreenMult_collect = chckPntIn.toScreenMult_collect;
    toScreenNoTreat_collect = chckPntIn.toScreenNoTreat_collect;
    toScreenNeg_collect = chckPntIn.toScreenNeg_collect; 
    toScreenTreat_collect = chckPntIn.toScreenTreat_collect; 
    toScreenTreatHystMult_collect = chckPntIn.toScreenTreatHystMult_collect; 
    
    % ART
    import java.util.LinkedList
    artDistList = chckPntIn.artDistList;
    artDist = chckPntIn.artDist;
    artTreatTracker = chckPntIn.artTreatTracker; %zeros(length(s) - 1 , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , gender , age , risk);
    artDiscont = chckPntIn.artDiscont;
end

%% Main body of simulation
% disp(['Simulating period from ' num2str(startYear) ' to ' num2str(endYear) ...
%     ' with ' num2str(stepsPerYear), ' steps per year.'])
disp(' ')
disp('Simulation running...')
disp(' ')
% progressbar('Simulation Progress')

for i = iStart : length(s) - 1
    year = startYear + s(i) - 1;
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
        fromNonHivAll = sort(toInd(allcomb(1 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1:risk))); 
        toHivAll = sort(toInd(allcomb(4 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1:risk)));
        
        % Distribute HIV infections 
        % % (HPV-)        
        % % popIn(fromNonHivNonHpv) = (1 - 0.002) .* popIn_init(fromNonHivNonHpv);    % reduce non-HIV infected
        % % popIn(toHivNonHpv) = (0.002) .* popIn_init(fromNonHivNonHpv);    % increase HIV infected ( male/female, age groups 4-6) (% prevalence)
        % % (HPV+)
        % % popIn(fromNonHivHpv) = (1 - 0.001) .* popIn_init(fromNonHivHpv);    % reduce non-HIV infected
        % % popIn(toHivHpv) = (0.001) .* popIn_init(fromNonHivHpv);    % increase HIV infected ( male/female, age groups 4-6) (% prevalence)
        popIn(fromNonHivAll) = (1 - 0.002) .* popIn_init(fromNonHivAll);    % reduce non-HIV infected
        popIn(toHivAll) = (0.002) .* popIn_init(fromNonHivAll);    % increase HIV infected ( male/female, age groups 4-6) (% prevalence)
        
    end

    if hpvOn % (Note: model not set up for HPV to be off; would need to add hpvOn bool logic to seeding of initial HPV infections)
        % HPV NATURAL HISTORY
        % Progression and clearance of HPV
        % Progression and regression of precancerous lesions
        % Development and progression of cervical cancer
        % Cervical cancer-associated mortality by stage and HIV status/CD4 count
        % Waning vaccine immunity
        [~ , pop , newCC(i , : , : , :) , ccDeath_untreat(i , : , : , :), ccDeath_treat(i , : , : , :), ...
            ccDeath_treat_stage(i , : , : , : , :)] ...
            = ode4xtra(@(t , pop) ...
            hpvCCNH(t , pop , hpv_hivClear , rImmuneHiv , c3c2Mults , c2c1Mults , c2c3Mults , c1c2Mults , muCC , muCC_ud , muCC_d , ...
            normalhpvVaxInds , immunehpvVaxInds , infhpvVaxInds , normalhpvNonVaxInds , ...
            immunehpvNonVaxInds , infhpvNonVaxInds , cin3hpvVaxIndsFrom , ccLochpvVaxIndsTo , ...
            ccLochpvVaxIndsFrom , ccReghpvVaxInds , ccDisthpvVaxInds , ...
            cin3hpvNonVaxIndsFrom , ccLochpvNonVaxIndsTo , ccLochpvNonVaxIndsFrom , ...
            ccReghpvNonVaxInds , ccDisthpvNonVaxInds , cin1hpvVaxInds , ...
            cin2hpvVaxInds , cin3hpvVaxInds , cin1hpvNonVaxInds , ...
            cin2hpvNonVaxInds , cin3hpvNonVaxInds , fromVaxNoScrnInds , ...
            fromVaxScrnInds , toNonVaxNoScrnInds , toNonVaxScrnInds , ...
            kInf_Cin1 , kCin1_Cin2 , kCin2_Cin3 , ...
            kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , kCin1_Inf , rNormal_Inf , ...
            rImmune , fImm , kRL , kDR , maleHpvClearMult , rVaxWane , disease , ...
            age , hpvVaxStates , hpvNonVaxStates , hpvTypeGroups , ccLochpvVaxIndsFrom_treat , ...
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
            % Symptomatic detection 

            [dPop , newScreen(i , : , : , : , : , : , :), ccTreat(i, : , : , : , : , : , : , :), ...
                ccSymp(i, : , : , : , : , : , : , :), ...
                toScreenMult_collect(i,:,:,:,:,:,:,:), toScreenNoTreat_collect(i,:,:,:,:,:,:,:), toScreenNeg_collect(i,:,:,:,:,:,:,:), ...
                toScreenTreat_collect(i,:,:,:,:,:,:,:), toScreenTreatHystMult_collect(i,:,:,:,:,:,:,:)] ...
                = hpvScreen(pop , ...
                    disease , viral , age , hpvVaxStates , hpvNonVaxStates , intervens, endpoints , risk , ...
                    screenYrs , screenAlgs , year , stepsPerYear , screenAgeAll , screenAgeS , ...
                    noVaxNoScreen , noVaxToScreen , vaxNoScreen , vaxToScreen , noVaxToScreenTreatImm , ...
                    vaxToScreenTreatImm , noVaxToScreenTreatImmVaxHpv , vaxToScreenTreatImmVaxHpv , ...
                    noVaxToScreenTreatImmNonVaxHpv , vaxToScreenTreatImmNonVaxHpv , ...        
                    noVaxToScreenTreatHpv , vaxToScreenTreatHpv , noVaxToScreenTreatVaxHpv , ...
                    vaxToScreenTreatVaxHpv , noVaxToScreenTreatNonVaxHpv , vaxToScreenTreatNonVaxHpv , ...
                    noVaxToScreenHyst , vaxToScreenHyst , numScreenAge , noVaxToScreenCancerNoTreat , noVaxToScreenCancerTreat , ...
                    vaxToScreenCancerNoTreat , vaxToScreenCancerTreat, hystMult, kSymp, udPop, udPopNoTreat, udPopTreat, udPopHyst, ...
                    vaxToScreenCancerNegScreen, noVaxToScreenCancerNegScreen);
            
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
            [dPop , ccSymp(i,:,:,:,:,:,:,:)] = symptomaticDetection(pop , ...
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
        newHpvNonVax(i , : , : , : , : , :) , newImmHpvNonVax(i , : , : , : , : , :) , newHiv(i , : , : , : , : , : , :)] = ...
        ode4xtra(@(t , pop) mixInfect(t , pop , ...
        stepsPerYear , year , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , intervens , gender , ...
        age , risk , fivYrAgeGrpsOn , hpvTypeGroups , ageSexDebut , gar , epsA_vec , epsR_vec , yr , ...
        partnersM , partnersF , ...
        beta_hpvVax_mod , beta_hpvNonVax_mod , vaxInds , nonVInds , ...
        lambdaMultImm , lambdaMultVax , artHpvMult , hpv_hivMult , ...
        hpvVaxSus , hpvVaxImm , hpvVaxInf , hpvNonVaxSus , hpvNonVaxImm , hpvNonVaxInf , ...
        circProtect , condProtect , condUse , betaHIV_mod , ...
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
        [~ , pop , hivDeaths(i , : , : , :) , artTreatTracker(i , : , : , : , : , :) , ...
            artDiscont(i , : , : , : , : , :)] =...
            ode4xtra(@(t , pop) hivNH(t , pop , vlAdvancer , muHIV , dMue , mue3 , mue4 , artDist , ...
            kCD4 ,  artYr_vec , artM_vec , artF_vec , minLim , maxLim , disease , viral , ...
            hpvVaxStates , hpvNonVaxStates , endpoints , gender , age , risk , ...
            ageSexDebut , hivInds , stepsPerYear , year) , tspan , popIn);
        popIn = pop(end , :);
        artTreat = artTreatTracker(i , : , : , : , : , :);
        artTreat = reshape(artTreat , [numel(artTreat) , 1]);
        artDistList.add(artTreat); %sum(sum(sum(artTreat , 3) , 4) , 5)
        if artDistList.size() >= stepsPerYear * 2
            artDistList.remove(); % remove CD4 and VL distribution info for people initiating ART more than 2 years ago
        end
        artDist = calcDist(artDistList , disease , viral , gender , age , risk);
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
            gender , age , fivYrAgeGrpsOn , fertMat , fertMat2 , fertMat3 , fertMat4 , ...
            hivFertPosBirth , hivFertNegBirth , hivFertPosBirth2 , hivFertNegBirth2 , ...
            hivFertPosBirth3 , hivFertNegBirth3 , hivFertPosBirth4 , hivFertNegBirth4 , ...
            dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ...
            dFertPos3 , dFertNeg3 , dFertMat3 , ...
            deathMat , deathMat2 , deathMat3 , deathMat4 , ...
            dDeathMat , dDeathMat2 , dDeathMat3 , ...
            MTCTRate , ageInd , riskInd , riskDist , ...
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
    
    if ((year >= vaxStartYear) && (vaxRate > 0))
        % HPV VACCINATION
        % School-based vaccination regimen
        [dPop , vaxdSchool(i , :)] = hpvVaxSchool(popIn , disease , viral , risk , ...
            hpvVaxStates , hpvNonVaxStates , endpoints , intervens , vaxG , vaxAge , ...
            vaxRate , toInd);
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
  
    if rem(year , 50) == 0.0
        savdir = [pwd , '/HHCoM_Results/'];
        save(fullfile(savdir , pathModifier, '') , 'fivYrAgeGrpsOn' , 'tVec' ,  'popVec' , 'newHiv' , ...
            'newHpvVax' , 'newImmHpvVax' , 'newHpvNonVax' , 'newImmHpvNonVax' , ...
            'hivDeaths' , 'deaths' , 'ccDeath' , 'menCirc' , 'vaxdSchool' , ...
            'newScreen' , 'ccDeath_treat', 'ccDeath_untreat', 'ccDeath_treat_stage', ...
            'newCC' , 'artDist' , 'artDistList' , 'artTreatTracker' , 'artDiscont' , ...
            'ccSymp' , 'ccTreat' , ...
            'startYear' , 'endYear' , 'i' , '-v7.3', ...
            'toScreenMult_collect', 'toScreenNoTreat_collect', 'toScreenNeg_collect', 'toScreenTreat_collect', 'toScreenTreatHystMult_collect' , 'kSymp');
    end

    disp(['Reached year ' num2str(year)])


end
popLast = sparse(popVec(end-1 , :));
disp(['Reached year ' num2str(endYear)])
popVec = sparse(popVec); % compress population vectors

%% Save results
savdir = [pwd , '/HHCoM_Results/'];
save(fullfile(savdir , pathModifier, '') , 'fivYrAgeGrpsOn' , 'tVec' ,  'popVec' , 'newHiv' , ...
    'newHpvVax' , 'newImmHpvVax' , 'newHpvNonVax' , 'newImmHpvNonVax' , ...
    'hivDeaths' , 'deaths' , 'ccDeath' , 'menCirc' , 'vaxdSchool' , ...
    'newScreen' , 'ccDeath_treat', 'ccDeath_untreat', 'ccDeath_treat_stage', ...
    'newCC' , 'artDist' , 'artDistList' , 'artTreatTracker' , 'artDiscont' , ...
    'ccSymp' , 'ccTreat' , ...
    'startYear' , 'endYear' , 'i' , 'popLast' , '-v7.3', ...
    'toScreenMult_collect', 'toScreenNoTreat_collect', 'toScreenNeg_collect', 'toScreenTreat_collect', 'toScreenTreatHystMult_collect' , 'kSymp');

disp(' ')
disp('Simulation complete.')
toc
% profile viewer

%% Calculate summed log-likelihood
% commented out temporarily for testing
% if calibBool    
%     negSumLogL = likeFun(popVec , newCC , cinPos2002_dObs , cinNeg2002_dObs ,...
%         hpv_hiv_dObs , hpv_hivNeg_dObs , hivPrevM_dObs , hivPrevF_dObs , ...
%         hpv_hivM2008_dObs , hpv_hivMNeg2008_dObs , ccInc2012_dObs , ccInc2018_dObs , cc_dist_dObs , ...
%         cin3_dist_dObs , cin1_dist_dObs , hpv_dist_dObs , ...
%         popAgeDist_dObs , totPopSize_dObs , cinPos2015_dObs , cinNeg2015_dObs , toInd , ...
%         disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , intervens , ...
%         age , gender , risk , startYear , stepsPerYear , annlz)
%     
%     if isnan(negSumLogL)
%         negSumLogL = -10000000.00;
%     end
%     
% else
%     negSumLogL = likeFun(popVec , newCC , cinPos2002_dObs , cinNeg2002_dObs ,...
%         hpv_hiv_dObs , hpv_hivNeg_dObs , hivPrevM_dObs , hivPrevF_dObs , ...
%         hpv_hivM2008_dObs , hpv_hivMNeg2008_dObs , ccInc2012_dObs , ccInc2018_dObs , cc_dist_dObs , ...
%         cin3_dist_dObs , cin1_dist_dObs , hpv_dist_dObs , ...
%         popAgeDist_dObs , totPopSize_dObs , cinPos2015_dObs , cinNeg2015_dObs , toInd , ...
%         disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , intervens , ...
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
