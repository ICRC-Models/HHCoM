% Future simulation module
% Accepts population vector from calibrated natural history model as input

function futureSim(calibBool , pIdx , paramsSub , paramSet , paramSetIdx , tstep_abc , date)    % input variables when using a calibration parameter set
% futureSim(0 , [] , [] , [] , [] , 0 , '19May20')    % input variables when running from command window using hand-calibrated, hard-coded parameter values
% Note: if you hard-code the "pathModifier" file output name variable below, then the date, paramSetIdx, and tstep_abc input values here are just dummy values and unused

% profile clear;

%% Cluster information
%pc = parcluster('local');    % create a local cluster object
%pc.JobStorageLocation = strcat('/gscratch/csde/carajb' , '/' , getenv('SLURM_JOB_ID'))    % explicitly set the JobStorageLocation to the temp directory that was created in the sbatch script
%parpool(pc , str2num(getenv('SLURM_CPUS_ON_NODE')))    % start the pool with max number workers

%%  Variables/parameters to set based on your scenario

% LOAD OUTPUT OF HISTORICAL SIMULATION AS INITIAL CONDITIONS FOR FUTURE SIMULATION
%historicalIn = load([pwd , '/HHCoM_Results/toNow_16Apr20_noBaseVax_baseScreen_hpvHIVcalib_0_1_test3_round1calib']); % ***SET ME***: name for historical run input file 
historicalIn = load([pwd , '/HHCoM_Results/toNow_' , date , '_baseVax057_baseScreen_baseVMMC_fertDec042-076_2020ARTfxd_trackCD4-Discont_discontFxd_DoART_S1_' , num2str(tstep_abc) , '_' , num2str(paramSetIdx)]); % ***SET ME***: name for historical run output file 

% DIRECTORY TO SAVE RESULTS
%pathModifier = '16Apr20_noBaseVax_baseScreen_hpvHIVcalib_0_1_test3_round1calib_050futureFert_WHOP1_SCES012'; % ***SET ME***: name for simulation output file
pathModifier = [date , '_baseVax057_baseScreen_baseVMMC_fertDec042-076-052_2020ARTfxd_trackCD4-Discont_discontFxd_diagHiv075_DoART_S2_' , num2str(tstep_abc) , '_' , num2str(paramSetIdx)]; % ***SET ME***: name for simulation output file
% Directory to save results
if ~ exist([pwd , '/HHCoM_Results/Vaccine' , pathModifier, '/'])
    mkdir ([pwd, '/HHCoM_Results/Vaccine' , pathModifier, '/'])
end

% AGE GROUPS
fivYrAgeGrpsOn = 1; % choose whether to use 5-year (fivYrAgeGrpsOn=1) or 1-year age groups (fivYrAgeGrpsOn=0)

% LAST YEAR
lastYear = 2061; % ***SET ME***: end year of simulation run

% SCREENING
% Instructions: Choose one screenAlgorithm, and modify the following screening parameters if appropriate.
screenAlgorithm = 1; % ***SET ME***: screening algorithm to use (1 for baseline, 2 for CISNET, 3 for WHOa, 4 for WHOb)
hivPosScreen = 0; % ***SET ME***: 0 applies same screening algorithm (screenAlgorithm) for all HIV states; 1 applies screenAlgorithm to HIV+ and screenAlgorithmNeg to HIV-
screenAlgorithmNeg = 1; % ***SET ME***: If hivPosScreen=1, screening algorithm to use for HIV- persons (1 for baseline, 2 for CISNET, 3 for WHOa, 4 for WHOb) 
whoScreenAges = [8 , 10]; %[6 , 7 , 8 , 9 , 10]; % ***SET ME***: ages that get screened when using the WHOa algorithm
whoScreenAgeMults = [0.20 , 0.20]; %[0.40 , 0.40 , 0.20 , 0.40 , 0.40]; % ***SET ME***: vector of equal length to whoScreenAges, fraction representing number of cohorts in each age range being screened

% VACCINATION
% Instructions: The model will run a scenario for each school-based vaccine coverage listed, plus a scenario with only baseline vaccine coverage.
%   If you want no vaccination in your baseline scenario, set baseline vaccine coverage to zero. The school-based vaccine coverage of each scenario is applied to all 
%   ages listed in that section. Therefore, if you assume baseline vaccination, your list of ages in the school-based vaccination algorithm should 
%   include the age of baseline vaccination, and school-based vaccine coverage should be at least baseline vaccine coverage.
%   If turned on, catch-up vaccine coverage is applied on top of all school-based vaccination scenarios, but not in the baseline vaccination only scenario. 
%   Distinct from the functionality of the school-based vaccination algorithm, catch-up vaccination coverage is defined by age group. Catch-up vaccination 
%   age groups should be exclusive of the school-based vaccination age groups.
%   If limited-vaccine years is turned on, this contraint is applied at the beginning of all the school-based vaccination scenarios, but not in the baseline 
%   vaccination only scenario. After the designated number of vaccine limited years has passed, the model will use the school based vaccination parameters 
%   and catch-up vaccination parameters if turned on.
% Example: 
%   Scenario 1: limited vaccine years --> school-based regimen for ages 9-14 at 86% coverage + catch-up coverage
%   Scenario 2: limited vaccine years --> school-based regimen for ages 9-14 at 90% coverage + catch-up coverage
%   Scenario 3: baseline regimen for age 9 at 86% coverage

% Common parameters
vaxEff = 1.0;  % 9v-vaccine efficacy, used for all vaccine regimens present
waning = 0;    % turn waning on or off

% Parameters for baseline vaccination regimen  % ***SET ME***: coverage for baseline vaccination of 9-year-old girls
vaxAgeB = [2];    % age groups to vaccinate
vaxCoverB = 0.57; %0.86;    % (9 year-old coverage * bivalent vaccine efficacy adjustment (0.7/0.9 proportion of cancers prevented) before 2020); last dose, first dose pilot
vaxGB = 2;   % indices of genders to vaccinate (1 or 2 or 1,2)

%Parameters for school-based vaccination regimen  % ***SET ME***: coverage for school-based vaccination of 9-14 year-old girls
vaxAge = [2 , 3];    % age groups to vaccinate
vaxCover = [0.57];    % vaccine coverages
vaxG = [2];   % indices of genders to vaccinate (1 or 2 or 1,2)

% Parameters for catch-up vaccination regimen
vaxCU = 0;    % turn catch-up vaccination on or off  % ***SET ME***: 0 for no catch-up vaccination, 1 for catch-up vaccination
hivPosVaxCU = 1;    % ***SET ME***: 0 applies catch-up vaccination algorithm for all HIV states; 1 applies catch-up vaccination only to HIV+ 
vaxAgeCU = [4 : 10];    % ages catch-up vaccinated % ***SET ME***: ages for catch-up vaccination
vaxCoverCU = [ones(1,length(vaxAgeCU)-1).*0.50 , 0.50*0.20];   % coverage for catch-up vaccination by ages catch-up vaccinated % ***SET ME***: coverage for catch-up vaccination by age, *adjustment factor if fraction of 5-year cohort
vaxGCU = [2];    % indices of genders to catch-up vaccinate (1 or 2 or 1,2)

% Parameters for vaccination during limited-vaccine years
vaxLimit = 0;    % turn vaccine limit on or off
vaxLimitYrs = 5;    % number years for which vaccines are limited
vaxLimitPerYr = 20000;    % total vaccines available per year for all interventions
vaxAgeL = 5;    % age group to vaccinate
vaxCoverL = 0.5;    % vaccine coverage
vaxGL = 2;    % index of gender to vaccinate during limited-vaccine years

% HIV TESTING CAMPAIGN
propHivDiagBaseline = [0.78 , 0.889]; % proportion diagnosed from SABSSMV (males, females)
propDiagOneYear = (1 - 0.41);
hivTestCampYrs = [2020 : 5 : lastYear-1];
hivTestCampCov = 0.75;

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
    c3c2Mults , c2c1Mults , muCC , kRL , kDR , artHpvMult , ...
    hpv_hivMult , maleHpvClearMult , ...
    condUse , screenYrs , hpvScreenStartYear , waning , ...
    artYr , maxRateM , maxRateF , ...
    artYr_vec , artM_vec , artF_vec , minLim , maxLim , ...
    circ_aVec , vmmcYr_vec , vmmc_vec , vmmcYr , vmmcRate , ...
    hivStartYear , circStartYear , circNatStartYear , vaxStartYear , ...
    baseline , cisnet , who , whob , circProtect , condProtect , MTCTRate , ...
    hyst , OMEGA , ...
    ccInc2012_dObs , cc_dist_dObs , cin3_dist_dObs , ...
    cin1_dist_dObs , hpv_dist_dObs , cinPos2002_dObs , cinNeg2002_dObs , ...
    hpv_hiv_dObs , hpv_hivNeg_dObs , hpv_hivM2008_dObs , hpv_hivMNeg2008_dObs , ...
    hivPrevM_dObs , hivPrevF_dObs , popAgeDist_dObs , totPopSize_dObs , ...
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
    infhpvNonVaxInds , ageInd , riskInd , deathInds , ...
    hivNegNonVMMCinds , hivNegVMMCinds , ...
    vlAdvancer , ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , fertMat3 , hivFertPosBirth3 , hivFertNegBirth3 , ...
    fertMat4 , hivFertPosBirth4 , hivFertNegBirth4 , ...
    dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ...
    dFertPos3 , dFertNeg3 , dFertMat3 , deathMat , deathMat2 , deathMat3 , deathMat4 , ...
    dDeathMat , dDeathMat2 , dDeathMat3 , dMue] = loadUp2(fivYrAgeGrpsOn , calibBool , pIdx , paramsSub , paramSet);

%% Screening

% WHO screening algorithm - version a
who.screenAge = whoScreenAges;
who.screenAgeMults = whoScreenAgeMults;

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
                            noVaxToScreenHyst(d,v,aS,r) = toInd(allcomb(d , v , 6 , 6 , 4 , 3 , 2 , a , r));
                            vaxToScreenHyst(d,v,aS,r) = toInd(allcomb(d , v , 6 , 6 , 4 , 4 , 2 , a , r));
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

%% Vaccination

% Set up differential HIV vaccination for catch-up vaccination regimen
if hivPosVaxCU
    vaxDiseaseIndsCU = [3 : 8];
else
    vaxDiseaseIndsCU = [1 : disease];
end

% Set up testParams vector for multiple school-based vaccination regimens
%   Example:
%   90% efficacy against 70% of CC types, 100% efficacy against 70% of types, 100% efficacy against 90% of types
%   vaxEff = [0.9 * 0.7 , 0.7 , 0.9]; 
testParams = allcomb(vaxCover , vaxEff); % test scenarios consist of all combinations of school-based vaccine coverage and efficacy
testParams = [testParams ; [vaxCoverB , vaxEff]]; % append baseline vaccination scenario to test scenarios
nTests = size(testParams , 1); % counts number of school-based scenarios to test + baseline scenario
testParams2(1:(nTests-1),1) = {vaxAge}; % age and gender to use with each school-based vaccination test scenario
testParams2(1:(nTests-1),2) = {vaxG};
testParams2(nTests,1) = {vaxAgeB}; % append baseline vaccination age and gender
testParams2(nTests,2) = {vaxGB};

if vaxCU
    vaxCoverCUmat = ones(nTests,length(vaxAgeCU)) .* vaxCoverCU;
    vaxCoverCUmat(end,:) = 0.0;
else
    vaxCoverCUmat = zeros(nTests,length(vaxAgeCU));    % have to declare these even if vaxCU=0 because parfor is dumb
end
if vaxLimit
    vaxCoverLmat = ones(nTests,1) .* vaxCoverL;
    vaxCoverLmat(end,1) = 0.0;
else
    vaxCoverLmat = zeros(nTests,1);    % have to declare these even if vaxLimit=0 because parfor is dumb
end

lambdaMultVaxMat = zeros(age , nTests); % age-based vector for modifying lambda based on vaccination status
vaxEffInd = repmat(1 : length(vaxEff) , 1 , (nTests) /length(vaxEff));
for n = 1 : nTests
    % No waning
    lambdaMultVaxMat(min(testParams2{n , 1}) : age , n) = vaxEff(vaxEffInd(n));
    
    % Waning
    effPeriod = 20; % number of years that initial efficacy level is retained
    wanePeriod = 20; % number of years over which initial efficacy level wanes
    if waning 
        % Following a period (in years) where original efficacy is retained, 
        % specified by 'effPeriod' , linearly scale down vaccine efficacy 
        % to 0% over time period specificed by 'wanePeriod'
        % To make waning rate equal in all scenarios, the linear rate of 
        % waning is based on the least effective initial vaccine efficacy scenario.        
        kWane = min(vaxEff) / round(wanePeriod / 5);     
        vaxInit = vaxEff(vaxEffInd(n));
        lambdaMultVaxMat(round(effPeriod / 5) + min(testParams2{n , 1}) - 1 : age , n) = ...
            max(0 , linspace(vaxInit , ...
            vaxInit - kWane * (1 + age - (round(wanePeriod / 5) + min(testParams2{n , 1}))) ,...
            age - (round(wanePeriod / 5) + min(testParams2{n , 1})) + 2)'); % ensures vaccine efficacy is >= 0
    end
end

%% Simulation
%profile on

for n = nTests
    simNum = n;
    vaxEff = testParams(n , 2);
    lambdaMultVax = 1 - lambdaMultVaxMat(: , n);
    vaxRate = testParams(n , 1);
    vaxAge = testParams2{n , 1};
    vaxG = testParams2{n , 2};
    if vaxCU
        vaxCoverCU = vaxCoverCUmat(n,:);
    end
    if vaxLimit
        vaxRemain = vaxLimitPerYr;
        vaxCoverL = vaxCoverLmat(n);
    end
    
    % Initial population
    popIn = historicalIn.popLast; % initial population to "seed" model
    
    % Initialize time vector
    yearsF = lastYear - currYear;
    s = 1 : timeStep : yearsF + 1;
    tVec = linspace(currYear , lastYear-timeStep , length(s)-1);
    
    % Initialize other vectors
    popVec = spalloc(length(s) - 1 , prod(dim) , 10 ^ 8);
    popVec(1 , :) = popIn;
    deaths = zeros(length(s) - 1 , disease , gender , age);
    newHiv = zeros(length(s) - 1 , hpvVaxStates , hpvNonVaxStates , endpoints , gender , age , risk);
    transCD4 = zeros(length(s) - 1 , disease , gender , age);
    hivDeaths = zeros(length(s) - 1 , disease , gender , age);
    % newHpvVax = zeros(length(s) - 1 , gender , disease , age , risk , intervens);
    % newImmHpvVax = newHpvVax;
    % newHpvNonVax = newHpvVax;
    % newImmHpvNonVax = newHpvVax;
    % newCC = zeros(length(s) - 1 , disease , age , hpvTypeGroups); % track by HPV type causal to CC
    % newCin1 = newCC;
    % newCin2 = newCC;
    % newCin3 = newCC;
    ccDeath = zeros(length(s) - 1 , disease , age , hpvTypeGroups);
    % newScreen = zeros(length(s) - 1 , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , numScreenAge , risk , 2);
    % newTreatImm = newScreen;
    % newTreatHpv = newScreen;
    % newTreatHyst = newScreen;
    menCirc = zeros(length(s) - 1 , 1);
    % vaxdLmtd = zeros(length(s) - 1 , 1);
    % vaxdSchool = vaxdLmtd;
    % vaxdCU = vaxdLmtd;
    
    % ART
    import java.util.LinkedList
    artDistList = historicalIn.artDistList;
    artDist = historicalIn.artDist;
    artTreatTracker = zeros(length(s) - 1 , disease , viral , gender , age , risk);
    artDiscont = zeros(length(s) - 1 , disease , viral , gender , age , risk);
    
    % Additional HIV testing campaigns
    nHivDiagVec = zeros(length(s) - 1 , gender);
    nHivUndiagVec = zeros(length(s) - 1 , gender);
    nTestedNeg = zeros(length(s) - 1 , gender);
    nTestedUndiag = zeros(length(s) - 1 , gender);
    propHivDiag = zeros(length(s) - 1 , gender);
    aged1519 = zeros(length(s) - 1 , gender);
    aged7579 = zeros(length(s) - 1 , gender);
    
    %% Main body of simulation
    for i = 2 : length(s) - 1
        year = currYear + s(i) - 1;
        tspan = [s(i) , s(i + 1)]; % evaluate diff eqs over one time interval
        popIn = popVec(i - 1 , :);
        
        nHivDiag = nHivDiagVec(i-1 , :);
        nHivUndiag = nHivUndiagVec(i-1 , :);
        
        % HIV testing initial conditions
        for g = 1 : gender
            nHivPos(1 , g) = sumall(popIn(hivInds(3 : 8 , 1 : viral , g , 4 : age , 1 : risk , :)));
            nHivNeg(1 , g) = sumall(popIn(hivInds(1 : 2 , 1 : viral , g , 4 : age , 1 : risk , :)));
            if i == 2
                nHivDiag(1 , g) = nHivPos(1 , g) * propHivDiagBaseline(1 , g);        
                nHivUndiag(1 , g) = nHivPos(1 , g) * (1 - propHivDiagBaseline(1 , g));
                propHivDiag(i , g) = propHivDiagBaseline(1 , g);
            end
        end
        %disp('TRACKED HIV VS. COMPARTMENT HIV AT START OF ITERATION')
        %(nHivDiag + nHivUndiag) - nHivPos
        
        % HIV testing, calculated rather than tracked in compartments
        if any(year == hivTestCampYrs) || (year == currYear + (1/stepsPerYear))
            [nTestedNeg(i , :) , nTestedUndiag(i , :) , propHivDiag(i , :) , nHivDiag , nHivUndiag] ...
                = hivTest(popIn , hivTestCampCov , nHivPos , nHivNeg , ...
                nHivUndiag , nHivDiag , hivInds , viral , gender , age , risk);
            %year
            propHivDiag(i , :);
            %%%num2str(propHivDiag(i , :))
        else
            propHivDiag(i , :) = nHivDiag(1 , :) ./ (nHivDiag(1 , :) + nHivUndiag(1 , :));
            %%%num2str(propHivDiag(i , :))
        end
        
        for g = 1 : gender
            nHivPos1(1 , g) = sumall(popIn(hivInds(3 : 8 , 1 : viral , g , 4 : age , 1 : risk , :)));
        end
        
        if hpvOn
            % Progression/regression from initial HPV infection to
            % precancer stages and cervical cancer. Differential CC
            % death by CC stage and HIV status/CD4 count.
            [~ , pop , ccDeath(i , : , : , :)] ...
                = ode4xtra(@(t , pop) ...
                hpvCCNH(t , pop , hpv_hivClear , rImmuneHiv , c3c2Mults , c2c1Mults , muCC , ...
                normalhpvVaxInds , immunehpvVaxInds , infhpvVaxInds , normalhpvNonVaxInds , ...
                immunehpvNonVaxInds , infhpvNonVaxInds , cin3hpvVaxIndsFrom , ccLochpvVaxIndsTo , ...
                ccLochpvVaxIndsFrom , ccReghpvVaxInds , ccDisthpvVaxInds , ...
                cin3hpvNonVaxIndsFrom , ccLochpvNonVaxIndsTo , ccLochpvNonVaxIndsFrom , ...
                ccReghpvNonVaxInds , ccDisthpvNonVaxInds , cin1hpvVaxInds , ...
                cin2hpvVaxInds , cin3hpvVaxInds , cin1hpvNonVaxInds , ...
                cin2hpvNonVaxInds , cin3hpvNonVaxInds , kInf_Cin1 , kCin1_Cin2 , kCin2_Cin3 , ...
                kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , kCin1_Inf , rNormal_Inf , ...
                rImmune , fImm , kRL , kDR , maleHpvClearMult , disease , age , hpvVaxStates , ...
                hpvNonVaxStates , hpvTypeGroups) , tspan , popIn);
            popIn = pop(end , :);  % for next module
            if any(pop(end , :) <  0)
                disp('After hpv')
                break
            end
            
            %disp('CC DEATH')
            for g = 1 : gender
                nHivPos2(1 , g) = sumall(popIn(hivInds(3 : 8 , 1 : viral , g , 4 : age , 1 : risk , :)));
            end
            %disp('Compartments vs. outputted change')
            %(nHivPos1(1,2)-nHivPos2(1,2)) - sumall(ccDeath(i , 3 : 8 , 4 : age , 1 : hpvTypeGroups))
                   
            if (year >= hpvScreenStartYear)
                [dPop] ...
                    = hpvScreen(popIn , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , risk , ...
                    screenYrs , screenAlgs , year , stepsPerYear , screenAgeAll , screenAgeS , ...
                    noVaxNoScreen , noVaxToScreen , vaxNoScreen , vaxToScreen , noVaxToScreenTreatImm , ...
                    vaxToScreenTreatImm , noVaxToScreenTreatHpv , vaxToScreenTreatHpv , ...
                    noVaxToScreenTreatVaxHpv , vaxToScreenTreatVaxHpv , noVaxToScreenTreatNonVaxHpv , ...
                    vaxToScreenTreatNonVaxHpv , noVaxToScreenHyst , vaxToScreenHyst , numScreenAge , ageMultsComb);
                pop(end , :) = pop(end , :) + dPop;
                popIn = pop(end , :);  % for next module
                if any(pop(end , :) <  0)
                    disp('After hpv screen')
                    break
                end
            end
        end
        
        % HIV and HPV mixing and infection module. Protective effects of condom
        % coverage, circumcision, ART, PrEP (not currently used) are accounted for. 
        [~ , pop , newHiv(i , : , : , : , : , : , :)] = ...
            ode4xtra(@(t , pop) mixInfect(t , pop , ...
            stepsPerYear , year , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , intervens , gender , ...
            age , risk , fivYrAgeGrpsOn , hpvTypeGroups , ageSexDebut , gar , epsA_vec , epsR_vec , yr , ...
            partnersM , partnersF , ...
            beta_hpvVax_mod , beta_hpvNonVax_mod , vaxInds , nonVInds , ...
            lambdaMultImm , lambdaMultVax , artHpvMult , hpv_hivMult , ...
            hpvVaxSus , hpvVaxImm , hpvVaxInf , hpvNonVaxSus , hpvNonVaxImm , hpvNonVaxInf , ...
            circProtect , condProtect , condUse , betaHIV_mod , ...
            hivSus , toHiv , hivCurr) , tspan , popIn);
        popIn = pop(end , :);
        if any(pop(end , :) < 0)
            disp('After mixInfect')
            break
        end
        
        %disp('NEW HIV')
        for g = 1 : gender
            nHivPos3(1 , g) = sumall(popIn(hivInds(3 : 8 , 1 : viral , g , 4 : age , 1 : risk , :)));
            nHivPosArt3(1 , g) = sumall(popIn(hivInds(8 , 1 : viral , g , 4 : age , 1 : risk , :)));
        end
        for a = 4 : age
            nHivPosArtAgeF(1 , a) = sumall(popIn(hivInds(8 , 1 : viral , 2 , a , 1 : risk , :)));
        end
        %disp('Compartments vs. outputed change')
        %[nHivPos3(1,:)-nHivPos2(1,:)] - [sumall(newHiv(i , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 , 4 : age , 1 : risk)) , sumall(newHiv(i , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 2 , 4 : age , 1 : risk))]
        
        % HIV module, CD4 Progression, VL progression, ART initiation/dropout,
        % excess HIV mortality
        if hivOn
            [~ , pop , hivDeaths(i , : , : , :) , artTreatTracker(i , : , : , : , : , :) , ...
                transCD4(i , : , : , :) , artDiscont(i , : , : , : , : , :)] =...
                ode4xtra(@(t , pop) hivNH(t , pop , propHivDiag(i , :) , vlAdvancer , muHIV , dMue , mue3 , mue4 , artDist , ... 
                kCD4 , artYr_vec , artM_vec , artF_vec , minLim , maxLim , disease , viral , ...
                hpvVaxStates , hpvNonVaxStates , endpoints , gender , age , risk , ...
                ageSexDebut , hivInds , stepsPerYear , year) , tspan , popIn);
                popIn = pop(end , :);    
            %artTreatTracker(i , : , : , : , : , :) = artTreat;
            artTreat = artTreatTracker(i , : , : , : , : , :);
            artTreat = reshape(artTreat , [numel(artTreat) , 1]);
            artDistList.add(artTreat); %sum(sum(sum(artTreat , 3) , 4) , 5)
            if artDistList.size() >= stepsPerYear * 2
                artDistList.remove(); % remove CD4 and VL distribution info for people initiating ART more than 2 years ago
            end
            artDist = calcDist(artDistList , disease , viral , gender , age , ...
                risk); % 2 year average CD4 and VL distribution at time of ART initiation. Details where ART dropouts return to.
            if any(pop(end , :) < 0)
                disp('After hiv')
                break
            end
        end
        
        for g = 1 : gender
            nHivPos4(1 , g) = sumall(popIn(hivInds(3 : 8 , 1 : viral , g , 4 : age , 1 : risk , :)));
            nHivPosArt4(1 , g) = sumall(popIn(hivInds(8 , 1 : viral , g , 4 : age , 1 : risk , :)));
        end
        %disp('All HIV: compartments vs. outputed change')
        %[nHivPos3(1,:)-nHivPos4(1,:)] - [sumall(hivDeaths(i , 1 : disease , 1 , 4 : age)) , sumall(hivDeaths(i , 1 : disease , 2 , 4 : age))]
        
        artChangeOut = [nHivPosArt4(1,:)-nHivPosArt3(1,:)] - [(sumall(artTreatTracker(i , 3 : disease , 1 : viral , 1 , 4 : age , 1 : risk))- sumall(artDiscont(i , 3 : disease , 1 : viral , 1 , 4 : age , 1 : risk)) - sumall(hivDeaths(i , 8 , 1 , 4 : age))) , ...
            (sumall(artTreatTracker(i , 3 : disease , 1 : viral , 2 , 4 : age , 1 : risk))- sumall(artDiscont(i , 3 : disease , 1 : viral , 2 , 4 : age , 1 : risk)) - sumall(hivDeaths(i , 8 , 2 , 4 : age)))];
        if abs(sumall(artChangeOut)) > 10
            disp('HIV DEATHS')
            disp('ART: compartments vs. outputted change')
            artChangeOut
            for a = 4 : age
                (sumall(popIn(hivInds(8 , 1 : viral , 2 , a , 1 : risk , :)))-nHivPosArtAgeF(1,a)) ...
                    - (sumall(artTreatTracker(i , 3 : disease , 1 : viral , 2 , a , 1 : risk)) - ...
                    sumall(artDiscont(i , 3 : disease , 1 : viral , 2 , a , 1 : risk)) - sumall(hivDeaths(i , 8 , 2 , a)))
            end
        end
        %disp('Female ART outputted change breakdown: ART init, ART discont, ART deaths')
        %sumall(artTreatTracker(i , 3 : disease , 1 : viral , 2 , 4 : age , 1 : risk))
        %sumall(artDiscont(i , 3 : disease , 1 : viral , 2 , 4 : age , 1 : risk))
        %sumall(hivDeaths(i , 8 , 2 , 4 : age))
        
        % Birth, aging, risk redistribution module
        [~ , pop , deaths(i , : , : , :) , aged1519(i , :) , aged7579(i , :)] = ode4xtra(@(t , pop) ...
            bornAgeDieRisk(t , pop , year , disease , viral , ...
            gender , age , risk , fivYrAgeGrpsOn , deathInds , fertMat , fertMat2 , fertMat3 , fertMat4 , ...
            hivFertPosBirth , hivFertNegBirth , hivFertPosBirth2 , hivFertNegBirth2 , ...
            hivFertPosBirth3 , hivFertNegBirth3 , hivFertPosBirth4 , hivFertNegBirth4 , ...
            dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ... 
            dFertPos3 , dFertNeg3 , dFertMat3 , ...
            deathMat , deathMat2 , deathMat3 , deathMat4 , ...
            dDeathMat , dDeathMat2 , dDeathMat3 , ...
            MTCTRate , hivInds , ageInd , riskInd , riskDist , ...
            stepsPerYear , currYear , agesComb , noVaxScreen , noVaxXscreen , ...
            vaxScreen , vaxXscreen , hpvScreenStartYear) , tspan , popIn);
        popIn = pop(end , :);
        if any(pop(end , :) < 0)
            disp('After bornAgeDieRisk')
            break
        end
        
        %disp('DEMOGRAPHICS')
        for g = 1 : gender
            nHivPos5(1 , g) = sumall(popIn(hivInds(3 : 8 , 1 : viral , g , 4 : age , 1 : risk , :)));
        end
        %disp('Compartments vs. outputted change')
        %[nHivPos5(1,:)-nHivPos4(1,:)] - [(aged1519(i , 1) - aged7579(i , 1) - sumall(deaths(i , hivInds(3 : 8 , 1 : viral , 1 , 4 : age , 1 : risk , :)))) , (aged1519(i , 2) - aged7579(i , 2) - sumall(deaths(i , hivInds(3 : 8 , 1 : viral , 2 , 4 : age , 1 : risk , :))))]
        
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
        
        % Estimate changes to proportion PLWHIV diagnosed/undiagnosed according to
        % population dynamics
        [nHivDiag , nHivUndiag] = hivTestPopStats(popIn , propDiagOneYear , ...
            nHivDiag , nHivUndiag , ...
            newHiv(i , : , : , : , : , : , :) , deaths(i , : , : , :) , ...
            hivDeaths(i , : , : , :) , ccDeath(i , : , : , :) , ...
            aged1519(i , :) , aged7579(i , :) , hivInds , disease , viral , ...
            hpvVaxStates , hpvNonVaxStates , endpoints , ...
            gender , age , risk , hpvTypeGroups);
        
        if (year >= vaxStartYear)
            % If within first vaxLimitYrs-many vaccine-limited years
            if vaxLimit && ((year - currYear) <= vaxLimitYrs)
                % HPV vaccination module- vaccine limited years
                [dPop , vaxRemain] = hpvVaxLmtd(popIn , year , vaxLimitPerYr , ...
                    disease , viral , risk , hpvVaxStates , hpvNonVaxStates , endpoints , ...
                    intervens , vaxCoverL , vaxRemain , vaxGL , toInd);
                pop(end , :) = pop(end , :) + dPop;
                popIn = pop(end , :);
                if any(pop(end , :) < 0)
                    disp('After hpvVaxLmtd')
                    break
                end
            
            % If vaccines are not limited
            else
                % HPV vaccination module- school-based vaccination regimen
                [dPop] = hpvVaxSchool(popIn , disease , viral , risk , ...
                    hpvVaxStates , hpvNonVaxStates , endpoints , intervens , vaxG , vaxAge , ...
                    vaxRate , toInd);
                pop(end , :) = pop(end , :) + dPop;
                popIn = pop(end , :);
                if any(pop(end , :) < 0)
                    disp('After hpvVaxSchool')
                    break
                end
                
                % If present, apply catch-up vaccination regimen
                if vaxCU
                    % HPV vaccination module- catch-up vaccination regimen
                    [dPop] = hpvVaxCU(popIn , viral , risk , ...
                        hpvVaxStates , hpvNonVaxStates , endpoints , intervens , vaxAgeCU , ...
                        vaxCoverCU , vaxGCU , vaxDiseaseIndsCU , toInd);
                    pop(end , :) = pop(end , :) + dPop;
                    if any(pop(end , :) < 0)
                        disp('After hpvVaxCU')
                        break
                    end
                end
            end
        end
       
        % add results to population vector
        popVec(i , :) = pop(end , :);
        
        nHivDiagVec(i , :) = nHivDiag;
        nHivUndiagVec(i , :) = nHivUndiag;
        
    end
    popLast = sparse(popVec(end , :));
    popVec = sparse(popVec); % compress population vectors
    
    filename = ['vaxSimResult' , num2str(simNum)];
    if waning
        filename = ['vaxWaneSimResult' , num2str(simNum)];
    end
    
    parsave(filename , fivYrAgeGrpsOn , tVec ,  popVec , newHiv , transCD4 , ...
        hivDeaths , deaths , ccDeath , ...
        menCirc , ...
        artDist , artDistList , artTreatTracker , artDiscont , ... 
        nHivDiagVec , nHivUndiagVec , nTestedNeg , nTestedUndiag , propHivDiag , ...
        currYear , lastYear , vaxRate , vaxEff , popLast , pathModifier);
end
disp('Done')

%profile viewer

%%
%vaxCEA(pathModifier)
