% Future simulation module
% Accepts population vector from calibrated natural history model as input

function futureSim(calibBool , pIdx , paramsSub , paramSet , paramSetIdx , tstep_abc , date , username)    % input variables when using a calibration parameter set
% futureSim(0 , [] , [] , [] , [] , 0 , '19May20' , 'carajb')    % input variables when running from command window using hand-calibrated, hard-coded parameter values
% Note: if you hard-code the "pathModifier" file output name variable below, then the date, paramSetIdx, and tstep_abc input values here are just dummy values and unused

% profile clear;

%% Cluster information (only use if running multiple vaccination scenarios in parallel. Cannot do this if using parfor loop in calib2_runMultFutSims.m to run multiple parameter sets in parallel)
%pc = parcluster('local');    % create a local cluster object
%pc.JobStorageLocation = strcat('/gscratch/csde/' , username , '/' , getenv('SLURM_JOB_ID'))    % explicitly set the JobStorageLocation to the temp directory that was created in the sbatch script
%parpool(pc , str2num(getenv('SLURM_CPUS_ON_NODE')))    % start the pool with max number workers

%%  Variables/parameters to set based on your scenario

% LOAD OUTPUT OF HISTORICAL SIMULATION AS INITIAL CONDITIONS FOR FUTURE SIMULATION
%historicalIn = load([pwd , '/HHCoM_Results/toNow_16Apr20_noBaseVax_baseScreen_hpvHIVcalib_0_1_test3_round1calib']);
historicalIn = load([pwd , '/HHCoM_Results/toNow_' , date , '_2v57BaseVax_spCytoScreen_shortName_noVMMChpv_discontFxd_screenCovFxd_hivInt2017_' , num2str(tstep_abc) , '_' , num2str(paramSetIdx)] , ...
    'popLast' , 'artDistList' , 'artDist'); % ***SET ME***: name for historical run output file 

% DIRECTORY TO SAVE RESULTS
%pathModifier = '16Apr20_noBaseVax_baseScreen_hpvHIVcalib_0_1_test3_round1calib_050futureFert_WHOP1_SCES012';
pathModifier = [date , '_2v57BaseVax_spCytoScreen_shortName_noVMMChpv_discontFxd_screenCovFxd_hivInt2017_SA-S0_' , num2str(tstep_abc) , '_' , num2str(paramSetIdx)]; % ***SET ME***: name for simulation output file
% Directory to save results
if ~ exist([pwd , '/HHCoM_Results/' , pathModifier, '/'])
    mkdir ([pwd, '/HHCoM_Results/' , pathModifier, '/'])
end

% AGE GROUPS
fivYrAgeGrpsOn = 1; % choose whether to use 5-year (fivYrAgeGrpsOn=1) or 1-year age groups (fivYrAgeGrpsOn=0)

% LAST YEAR
lastYear = 2122; % ***SET ME***: end year of simulation run

% SCREENING
% Instructions: Choose one screenAlgorithm, and modify the following screening parameters if appropriate.
%   For example, if you want persons across all HIV states to follow the same screening pattern,
%   use sceScreenHivGrps={[1:8]} sceScreenAges={[8 , 10]} for 2x screening regardless of HIV status. 
%   If you want screening pattern to differ by HIV status, use sceScreenHivGrps={[1 : 2] , [3 : 8]} 
%   to designate different patterns for HIV-negative and HIV-positive women and 
%   sceScreenAges={[8 , 10] , [6 , 7 , 8 , 9 , 10]} for 2x screening among HIV-negative women and screening 
%   every 3 years among HIV-positive women.
screenAlgorithm = 3; % ***SET ME***: screening algorithm to use (1 for baseline, 2 for WHO, 3 for spCyto, 4 for spHpvDna, 5 for spGentyp, 6 for spAve , 7 for spHpvAve)
sceScreenCover = [0.0; 0.18; 0.48; 0.48;     0.48; 0.48; 0.48]; % Coverage over time (Years: [2000; 2003; 2016; currYear;     2023; 2030; 2045])
sceScreenHivGrps = {[1 : 8]}; % ***SET ME***: Groupings of HIV states with different screening ages
sceScreenAges = {[8]}; % ***SET ME***: screening ages that correspond to HIV state groupings

% VACCINATION
% Instructions: The model will set up a scenario for each school-based vaccine coverage listed in "vaxCover", plus a scenario with only baseline vaccine coverage as in "vaxCoverB".
%   If you want no vaccination in your baseline scenario, set baseline vaccine coverage to zero. The school-based vaccine coverage of each scenario is applied to all 
%   ages listed in that section. Therefore, if you assume baseline vaccination, your list of ages in the school-based vaccination algorithm should 
%   include the age of baseline vaccination, and school-based vaccine coverage should be at least baseline vaccine coverage.
%   If turned on, catch-up vaccine coverage is applied on top of all school-based vaccination scenarios, but not in the baseline vaccination only scenario. 
%   Distinct from the functionality of the school-based vaccination algorithm, catch-up vaccination coverage is defined by age group. Catch-up vaccination 
%   age groups should be exclusive of the school-based vaccination age groups.
%   If limited-vaccine years is turned on, this contraint is applied at the beginning of all the school-based vaccination scenarios, but not in the baseline 
%   vaccination only scenario. After the designated number of vaccine limited years has passed, the model will use the school based vaccination parameters 
%   and catch-up vaccination parameters if turned on.
% Example (scenarios set up for vaxCover = [0.86, 0.90]; and vaxCoverB = 0.86): 
%   Scenario 1: limited vaccine years --> school-based regimen for ages 9-14 at 86% coverage + catch-up coverage (to run scenario, set vaxCoverInd = 1)
%   Scenario 2: limited vaccine years --> school-based regimen for ages 9-14 at 90% coverage + catch-up coverage (to run scenario, set vaxCoverInd = 2)
%   Scenario 3: baseline regimen for age 9 at 86% coverage (to run scenario, set vaxCoverInd = 3)

% Common parameters
vaxEff = 1.0;  % 9v-vaccine efficacy, used for all vaccine regimens present
rVaxWane = 0.0; % rate of waning vaccine immunity

% Parameters for baseline vaccination regimen  % ***SET ME***: coverage for baseline vaccination of 9-year-old girls
vaxAgeB = [2];    % age groups to vaccinate
vaxCoverB = 0.0; %0.57; %0.86;    % (9 year-old coverage * bivalent vaccine efficacy adjustment (0.7/0.9 proportion of cancers prevented) before 2020); last dose, first dose pilot
vaxGB = 2;   % indices of genders to vaccinate (1 or 2 or 1,2); set stepsPerYear=8 in loadUp2.m if including vaccination of boys 

%Parameters for school-based vaccination regimen  % ***SET ME***: coverage for school-based vaccination of 9-14 year-old girls
vaxAge = [2 , 3];    % age groups to vaccinate
vaxCover = [0.57];    % vaccine coverages
vaxCoverInd = 1;    % index for the coverage in vaxCover vec to use for this simulation; use length(vaxCover)+1 to run the baseline scenario (Ex: vaxCoverInd=1 for specified scenario, vaxCoverInd=2 for baseline scenario)
vaxG = [2];   % indices of genders to vaccinate (1 or 2 or 1,2); set stepsPerYear=8 in loadUp2.m if including vaccination of boys 

% Parameters for catch-up vaccination regimen
vaxCU = 0;    % turn catch-up vaccination on or off  % ***SET ME***: 0 for no catch-up vaccination, 1 for catch-up vaccination
hivPosVaxCU = 1;    % ***SET ME***: 0 applies catch-up vaccination algorithm for all HIV states; 1 applies catch-up vaccination only to HIV+ 
vaxAgeCU = [4 : 5];    % ages catch-up vaccinated % ***SET ME***: ages for catch-up vaccination
vaxCoverCU = [ones(1,length(vaxAgeCU)).*0.90];   % coverage for catch-up vaccination by ages catch-up vaccinated % ***SET ME***: coverage for catch-up vaccination by age, *adjustment factor if fraction of 5-year cohort
vaxGCU = [2];    % indices of genders to catch-up vaccinate (1 or 2 or 1,2)

% Parameters for vaccination during limited-vaccine years
vaxLimit = 0;    % turn vaccine limit on or off
vaxLimitYrs = 5;    % number years for which vaccines are limited
vaxLimitPerYr = 20000;    % total vaccines available per year for all interventions
vaxAgeL = 5;    % age group to vaccinate
vaxCoverL = 0.5;    % vaccine coverage
vaxGL = 2;    % index of gender to vaccinate during limited-vaccine years

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
    c3c2Mults , c2c1Mults , c2c3Mults , c1c2Mults , muCC , kRL , kDR , artHpvMult , ...
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
    hivNegNonVMMCinds , hivNegVMMCinds , ...
    vlAdvancer , ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , fertMat3 , hivFertPosBirth3 , hivFertNegBirth3 , ...
    fertMat4 , hivFertPosBirth4 , hivFertNegBirth4 , ...
    dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ...
    dFertPos3 , dFertNeg3 , dFertMat3 , deathMat , deathMat2 , deathMat3 , deathMat4 , ...
    dDeathMat , dDeathMat2 , dDeathMat3 , dMue] = loadUp2(fivYrAgeGrpsOn , calibBool , pIdx , paramsSub , paramSet);

%% Screening
if (screenAlgorithm == 1)
    % Baseline screening algorithm
    screenAlgs = baseline;
    screenAlgs.genTypBool = 0;
elseif (screenAlgorithm == 2)
    % WHO screening algorithm
    screenAlgs = who;
    screenAlgs.genTypBool = 0;
elseif (screenAlgorithm == 3)
    % Screening paper cytology algorithm
    screenAlgs = spCyto;
    screenAlgs.genTypBool = 0;
elseif (screenAlgorithm == 4)
    % Screening paper HPV DNA -and-treat algorithm
    screenAlgs = spHpvDna;
    screenAlgs.genTypBool = 0;
elseif (screenAlgorithm == 5)
    % Screening paper HPV DNA+genotyping -and-treat algorithm
    screenAlgs = spGentyp;
    screenAlgs.genTypBool = 1;
elseif (screenAlgorithm == 6)
    % Screening paper AVE -and-treat algorithm
    screenAlgs = spAve;
    screenAlgs.genTypBool = 0;
elseif (screenAlgorithm == 7)
    % Screening paper HPV DNA + AVE triage -and-treat algorithm
    screenAlgs = spHpvAve;
    screenAlgs.genTypBool = 0;
end
screenAlgs.screenHivGrps = sceScreenHivGrps;
screenAlgs.screenAge = sceScreenAges;
screenAlgs.screenCover = sceScreenCover;
screenAlgs.screenCover_vec = cell(size(screenYrs , 1) - 1, 1); % save data over time interval in a cell array
for i = 1 : size(screenYrs , 1) - 1          % interpolate dnaTestCover values at steps within period
    period = [screenYrs(i) , screenYrs(i + 1)];
    screenAlgs.screenCover_vec{i} = interp1(period , screenAlgs.screenCover(i : i + 1 , 1) , ...
        screenYrs(i) : timeStep : screenYrs(i + 1));
end

% Create screening indices
numScreenAge = 0;
agesComb = [];
for n = 1 : length(sceScreenAges)
    numScreenAge = numScreenAge + length(sceScreenAges{n});
    agesComb = [agesComb , sceScreenAges{n}]; 
end
screenAgeAll = zeros(disease , viral , numScreenAge , risk , hpvVaxStates*hpvNonVaxStates*endpoints*intervens);
screenAgeS = zeros(disease , viral , numScreenAge , risk , hpvVaxStates*hpvNonVaxStates*endpoints*2);
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

for aS = 1 : numScreenAge
    a = agesComb(aS);
    
    for d = 1 : disease
        for v = 1 : viral
            for h = 1 : hpvVaxStates
                for s = 1 : hpvNonVaxStates
                    for x = 1 : endpoints
                        for r = 1 : risk
                            screenAgeAll(d,v,aS,r,:) = toInd(allcomb(d , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 2 , a , r)); 
                            screenAgeS(d,v,aS,r,:) = toInd(allcomb(d , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 3 : intervens , 2 , a , r));

                            noVaxNoScreen(d,v,h,s,x,aS,r) = sort(toInd(allcomb(d , v , h , s , x , 1 , 2 , a , r)));
                            noVaxToScreen(d,v,h,s,x,aS,r) = sort(toInd(allcomb(d , v , h , s , x , 3 , 2 , a , r)));
                            vaxNoScreen(d,v,h,s,x,aS,r) = sort(toInd(allcomb(d , v , h , s , x , 2 , 2 , a , r)));
                            vaxToScreen(d,v,h,s,x,aS,r) = sort(toInd(allcomb(d , v , h , s , x , 4 , 2 , a , r)));

                            noVaxToScreenTreatImm(d,v,aS,r) = toInd(allcomb(d , v , 7 , 7 , 1 , 3 , 2 , a , r));
                            vaxToScreenTreatImm(d,v,aS,r) = toInd(allcomb(d , v , 7 , 7 , 1 , 4 , 2 , a , r));
                            noVaxToScreenTreatImmVaxHpv(d,v,s,aS,r) = toInd(allcomb(d , v , 7 , s , 1 , 3 , 2 , a , r));
                            vaxToScreenTreatImmVaxHpv(d,v,s,aS,r) = toInd(allcomb(d , v , 7 , s , 1 , 4 , 2 , a , r));
                            noVaxToScreenTreatImmNonVaxHpv(d,v,h,aS,r) = toInd(allcomb(d , v , h , 7 , 1 , 3 , 2 , a , r));
                            vaxToScreenTreatImmNonVaxHpv(d,v,h,aS,r) = toInd(allcomb(d , v , h , 7 , 1 , 4 , 2 , a , r));
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
    lambdaMultVaxMat(min(testParams2{n , 1}) : age , n) = vaxEff(vaxEffInd(n));
end

%% Simulation
%profile on
n = vaxCoverInd; %parfor n = 1 : nTests (can only use parfor loop if not running multiple parameter sets in parallel
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
    deaths = zeros(length(s) - 1 , 1); %zeros(size(popVec));
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
    newScreen = zeros(length(s) - 1 , disease , hpvVaxStates , hpvNonVaxStates , endpoints , numScreenAge , 2);
    % newTreatImm = newScreen;
    % newTreatHpv = newScreen;
    % newTreatHyst = newScreen;
    menCirc = zeros(length(s) - 1 , 1);
    vaxdLmtd = zeros(length(s) - 1 , 1);
    vaxdSchool = vaxdLmtd;
    vaxdCU = vaxdLmtd;
    % ART
    import java.util.LinkedList
    artDistList = historicalIn.artDistList;
    artDist = historicalIn.artDist;
    artTreatTracker = zeros(length(s) - 1 , disease , viral , gender , age , risk);
    artDiscont = zeros(length(s) - 1 , disease , viral , gender , age , risk);
    
    %% Main body of simulation
    for i = 2 : length(s) - 1
        year = currYear + s(i) - 1;
        tspan = [s(i) , s(i + 1)]; % evaluate diff eqs over one time interval
        popIn = popVec(i - 1 , :);
        
        if hpvOn
            % HPV NATURAL HISTORY
            % Progression and clearance of HPV
            % Progression and regression of precancerous lesions
            % Development and progression of cervical cancer
            % Cervical cancer-associated mortality by stage and HIV status/CD4 count
            % Waning vaccine immunity
            [~ , pop , newCC(i , : , : , :) , ccDeath(i , : , : , :)] ...
                = ode4xtra(@(t , pop) ...
                hpvCCNH(t , pop , hpv_hivClear , rImmuneHiv , c3c2Mults , c2c1Mults , c2c3Mults , c1c2Mults , muCC , ...
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
                age , hpvVaxStates , hpvNonVaxStates , hpvTypeGroups) , tspan , popIn);
            popIn = pop(end , :);  % for next module
            if any(pop(end , :) <  0)
                disp('After hpv')
                break
            end
            
            if (year >= hpvScreenStartYear)
                % CERVICAL CANCER SCREENING AND TREATMENT
                % Screening
                % Treatment
                [dPop , newScreen(i , : , : , : , : , : , :)] ...
                    = hpvScreen(popIn , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , risk , ...
                    screenYrs , screenAlgs , year , stepsPerYear , screenAgeAll , screenAgeS , ...
                    noVaxNoScreen , noVaxToScreen , vaxNoScreen , vaxToScreen , noVaxToScreenTreatImm , ...
                    vaxToScreenTreatImm , noVaxToScreenTreatImmVaxHpv , vaxToScreenTreatImmVaxHpv , ...
                    noVaxToScreenTreatImmNonVaxHpv , vaxToScreenTreatImmNonVaxHpv , ...        
                    noVaxToScreenTreatHpv , vaxToScreenTreatHpv , noVaxToScreenTreatVaxHpv , ...
                    vaxToScreenTreatVaxHpv , noVaxToScreenTreatNonVaxHpv , vaxToScreenTreatNonVaxHpv , ...
                    noVaxToScreenHyst , vaxToScreenHyst , numScreenAge);
                pop(end , :) = pop(end , :) + dPop;
                popIn = pop(end , :);  % for next module
                if any(pop(end , :) <  0)
                    disp('After hpv screen')
                    break
                end
            end
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
        popIn = pop(end , :);
        if any(pop(end , :) < 0)
            disp('After mixInfect')
            break
        end
        
        % HIV NATURAL HISTORY
        % CD4 progression
        % Viral load progression
        % ART initiation, dicontinuation, and scale-up by CD4 count
        % HIV-associated mortality
        if hivOn
            [~ , pop , hivDeaths(i , : , : , :) , artTreatTracker(i , : , : , : , : , :) , ...
                artDiscont(i , : , : , : , : , :)] =...
                ode4xtra(@(t , pop) hivNH(t , pop , vlAdvancer , muHIV , dMue , mue3 , mue4 , artDist , ... 
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
        
        % HPV VACCINATION
        if (year >= vaxStartYear)
            % If within first vaxLimitYrs-many vaccine-limited years
            if vaxLimit && ((year - currYear) <= vaxLimitYrs)
                % HPV vaccination module- vaccine limited years
                [dPop , vaxdLmtd(i , :) , vaxRemain] = hpvVaxLmtd(popIn , year , vaxLimitPerYr , ...
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
                % School-based vaccination regimen
                [dPop , vaxdSchool(i , :)] = hpvVaxSchool(popIn , disease , viral , risk , ...
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
                    [dPop , vaxdCU(i , :)] = hpvVaxCU(popIn , viral , risk , ...
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
    end
    popLast = sparse(popVec(end , :));
    popVec = sparse(popVec); % compress population vectors

    filename = ['vaxSimResult' , num2str(simNum)];
    parsave(filename , fivYrAgeGrpsOn , tVec ,  popVec , newHiv ,...
        newHpvVax , newImmHpvVax , newHpvNonVax , newImmHpvNonVax , ...
        hivDeaths , deaths , ccDeath , ...
        newCC , menCirc , vaxdLmtd , vaxdSchool , vaxdCU , newScreen , ...
        artDist , artDistList , artTreatTracker , artDiscont , ... 
        currYear , lastYear , vaxRate , vaxEff , popLast , pathModifier);
%end
disp('Done')

%profile viewer

%%
%vaxCEA(pathModifier)
