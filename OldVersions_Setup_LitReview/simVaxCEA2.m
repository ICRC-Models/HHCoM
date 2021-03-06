% Future simulation module
% Accepts population vector from calibrated natural history model as input

%%
close all; clear all; clc
% profile clear;

%% Cluster information
% pc = parcluster('local');    % create a local cluster object
% pc.JobStorageLocation = strcat('/gscratch/csde/carajb' , '/' , getenv('SLURM_JOB_ID'))    % explicitly set the JobStorageLocation to the temp directory that was created in the sbatch script
% parpool(pc , str2num(getenv('SLURM_CPUS_ON_NODE')))    % start the pool with max number workers

%%  Variables/parameters to set based on your scenario

% LOAD POPULATION
popIn = load([pwd , '/HHCoM_Results/toNow_101819_5yrAgeGrps_noBaseVax_baseScreen_nonVhpv']); % ***SET ME***: name for historical run input file 
currPop = popIn.popLast;
artDistIn = popIn.artDist;
artDistListIn = popIn.artDistList;

% DIRECTORY TO SAVE RESULTS
pathModifier = '101819_5yrAgeGrps_noBaseVax_baseScreen_nonVhpv'; % ***SET ME***: name for simulation output file
if ~ exist([pwd , '/HHCoM_Results/Vaccine' , pathModifier, '/'])
    mkdir ([pwd, '/HHCoM_Results/Vaccine' , pathModifier, '/'])
end

% AGE GROUPS
fivYrAgeGrpsOn = 1; % choose whether to use 5-year or 1-year age groups

% LAST YEAR
lastYear = 2121; % ***SET ME***: end year of simulation run

% SCREENING
screenAlgorithm = 1; % ***SET ME***: screening algorithm to use (1 for baseline, 2 for CISNET, 3 for WHOa, 4 for WHOb)
hivPosScreen = 0; % ***SET ME***: 0 applies same screening algorithm (screenAlgorithm) for all HIV states; 1 applies screenAlgorithm to HIV+ and screenAlgorithmNeg to HIV-
screenAlgorithmNeg = 1; % ***SET ME***: If hivPosScreen=1, screening algorithm to use for HIV- persons (1 for baseline, 2 for CISNET, 3 for WHOa, 4 for WHOb) 
whoScreenAges = [8 , 10]; %[6 , 7 , 8 , 9 , 10]; %[26 , 29 , 32 , 35 , 38 , 41 , 44 , 47 , 50]; % ***SET ME***: ages that get screened when using the WHOa algorithm
whoScreenAgeMults = [0.20 , 0.20]; %[0.40 , 0.40 , 0.20 , 0.40 , 0.40];

% VACCINATION
vaxEff = [0.9];    % 9v-vaccine, used for all vaccine regimens present
waning = 0;    % turn waning on or off

% Parameters for baseline vaccination regimen  % ***SET ME***: coverage for baseline vaccination of 9-year-old girls
vaxAgeB = [2];
vaxCoverB = 0.0; %0.86*(0.7/0.9);    % (9 year-old coverage * bivalent vaccine efficacy adjustment)
vaxGB = 2;   % indices of genders to vaccinate (1 or 2 or 1,2)

%Parameters for school-based vaccination regimen  % ***SET ME***: coverage for school-based vaccination of 9-14 year-old girls
vaxAge = [2 , 3];
vaxCover = [0.8 , 0.9];
vaxG = [2];   % indices of genders to vaccinate (1 or 2 or 1,2)

% Parameters for catch-up vaccination regimen
vaxCU = 0;    % turn catch-up vaccination on or off  % ***SET ME***: 0 for no catch-up vaccination, 1 for catch-up vaccination
hivPosVaxCU = 1; % ***SET ME***: 0 applies catch-up vaccination algorithm for all HIV states; 1 applies catch-up vaccination only to HIV+ 
vaxAgeCU = [4 : 10]; %[16 : 27];    % ages catch-up vaccinated % ***SET ME***: ages for catch-up vaccination
vaxCoverCU = [ones(1,length(vaxAgeCU)-1).*0.50 , 0.50*0.20]; %0.50 % coverage for catch-up vaccination by ages catch-up vaccinated % ***SET ME***: coverage for catch-up vaccination by age
vaxGCU = [2];    % indices of genders to catch-up vaccinate (1 or 2 or 1,2)

% Parameters for vaccination during limited-vaccine years
vaxLimit = 0;    % turn vaccine limit on or off
vaxLimitYrs = 5;    % number years for which vaccines are limited
vaxLimitPerYr = 20000;    % total vaccines available per year for all interventions
vaxAgeL = 5;
vaxCoverL = 0.5;
vaxGL = 2;    % index of gender to vaccinate during limited-vaccine years

%% Model specs
% choose whether to model background hysterectomy
hyst = 0; % NOT UPDATED!!!!!!!!!!!!!!!!!
% choose whether to model HIV
hivOn = 1;
% choose whether to model HPV
hpvOn = 1;
if hpvOn
    disp('HPV module activated')
end
if hivOn
    disp('HIV module activated')
end

%% Save pre-loaded parameters and pre-calculated indices and matrices
loadUp2(fivYrAgeGrpsOn)

%% Load saved parameters
disp('Initializing. Standby...')

paramDir = [pwd , '/Params/'];

% Load saved parameters
load([paramDir, 'generalParams'], 'stepsPerYear' , 'timeStep' , ...
    'disease' , 'viral' , 'hpvVaxStates' , 'hpvNonVaxStates' , 'endpoints' , ...
    'intervens' , 'gender' , 'age' , 'risk' , 'hpvTypeGroups' , 'dim' , 'k' , 'toInd' , ...
    'sumall');
load([paramDir, 'demoBehavParams'], 'ageSexDebut' , 'mInit' , 'fInit' , 'partnersM' , 'partnersF' , ...
    'maleActs' , 'femaleActs' , 'riskDist' , 'mue' , 'epsA_vec' , 'epsR_vec' , 'yr');
load([paramDir, 'hivParams'], 'betaHIVM2F' , 'betaHIVF2M' , 'muHIV' , 'kVl' , 'kCD4');
load([paramDir, 'hpvParams'], 'perPartnerHpv_vax' , 'perPartnerHpv_nonV' , ...
    'fImm' , 'rImmune' , 'kCin1_Inf' , 'kCin2_Cin1' , 'kCin3_Cin2' , 'kCC_Cin3' , ...
    'rNormal_Inf' , 'kInf_Cin1' , 'kCin1_Cin2' , 'kCin2_Cin3' , 'lambdaMultImm' , ...
    'hpv_hivClear' , 'rImmuneHiv' , 'c3c2Mults' , 'c2c1Mults' , 'muCC' , ...
    'kRL' , 'kDR' , 'artHpvMult' , 'hpv_hivMult');
load([paramDir, 'intervenParams'], 'circ' , 'condUse' , ...
    'maxRateM1' , 'maxRateF1' , 'maxRateM2' , 'maxRateF2' , 'hivStartYear' , ...
    'circStartYear' , 'vaxStartYear' , 'baseline' , 'cisnet' , 'who' , 'whob' , ...
    'circProtect' , 'condProtect' , 'MTCTRate');
load([paramDir , 'calibData'], 'cinPos2002_obs' , 'cinNeg2002_obs' , ...
    'hpv_hiv_obs' , 'hpv_hivNeg_obs' , 'hpv_hivM2008_obs' , 'hpv_hivMNeg2008_obs' , ...
    'hivPrevM_obs' , 'hivPrevF_obs');

%% Load saved indices
load([paramDir,'mixInfectIndices'], 'mCurr' , 'fCurr' , 'mCurrArt' , 'fCurrArt' , ...
    'gar' , 'hivSus' , 'hpvVaxSus' , 'hpvVaxImm' , ...
    'hpvNonVaxSus' , 'hpvNonVaxImm' , ...
    'toHiv' , 'vaxInds' , 'nonVInds' , 'hpvVaxInf' , 'hpvNonVaxInf');
load([paramDir,'hivIndices'], 'hivInds');
load([paramDir,'hpvIndices'], 'cin3hpvVaxIndsFrom' , 'ccLochpvVaxIndsTo' , ...
    'ccLochpvVaxIndsFrom' , 'ccReghpvVaxInds' , 'ccDisthpvVaxInds' , ...
    'cin3hpvNonVaxIndsFrom' , 'ccLochpvNonVaxIndsTo' , 'ccLochpvNonVaxIndsFrom' , ...
    'ccReghpvNonVaxInds' , 'ccDisthpvNonVaxInds' , 'cin1hpvVaxInds' , ...
    'cin2hpvVaxInds' , 'cin3hpvVaxInds' , 'cin1hpvNonVaxInds' , ...
    'cin2hpvNonVaxInds' , 'cin3hpvNonVaxInds' , 'normalhpvVaxInds' , 'immunehpvVaxInds' , ...
    'infhpvVaxInds' , 'normalhpvNonVaxInds' , 'immunehpvNonVaxInds' , 'infhpvNonVaxInds');
load([paramDir,'ageRiskInds'], 'ageInd' , 'riskInd');

%% Load saved matrices
load([paramDir,'vlAdvancer'], 'vlAdvancer')
load([paramDir,'fertMat'], 'fertMat')
load([paramDir,'hivFertMats'], 'hivFertPosBirth' , 'hivFertNegBirth')
load([paramDir,'fertMat2'], 'fertMat2')
load([paramDir,'hivFertMats2'], 'hivFertPosBirth2' , 'hivFertNegBirth2')
load([paramDir,'deathMat'], 'deathMat')
load([paramDir,'circMat'], 'circMat')
load([paramDir,'circMat2'] , 'circMat2')

%% Time
c = fix(clock);
currYear = 2020; % c(1); % get the current year
years = lastYear - currYear;

%% Screening
screenYrs = [2000; 2003; 2016; currYear; 2023; 2030; 2045];
hpvScreenStartYear = screenYrs(1);

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
testParams = [testParams ; [vaxCoverB , vaxEff]]; % Append baseline vaccination scenario to test scenarios
nTests = size(testParams , 1); % counts number of school-based scenarios to test
testParams2(1:(nTests-1),1) = {vaxAge};
testParams2(1:(nTests-1),2) = {vaxG};
testParams2(nTests,1) = {vaxAgeB};
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

parfor n = 1 : nTests
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
    
    % Initialize time vector
    s = 1 : timeStep : years + 1;
    tVec = linspace(currYear , lastYear-timeStep , length(s)-1);
    
    % Initialize other vectors
    popVec = spalloc(length(s) - 1 , prod(dim) , 10 ^ 8);
    popIn = currPop; % initial population to "seed" model
    popVec(1 , :) = popIn;
    deaths = zeros(size(popVec));
    newHiv = zeros(length(s) - 1 , gender , age , risk);
    hivDeaths = zeros(length(s) - 1 , gender , age);
    newHpvVax = zeros(length(s) - 1 , gender , disease , age , risk , intervens);
    newImmHpvVax = newHpvVax;
    newHpvNonVax = newHpvVax;
    newImmHpvNonVax = newHpvVax;
    newCC = zeros(length(s) - 1 , disease , age , hpvTypeGroups); % track by HPV type causal to CC
%     newCin1 = newCC;
%     newCin2 = newCC;
%     newCin3 = newCC;
    ccDeath = newCC;
    newScreen = zeros(length(s) - 1 , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , numScreenAge , risk , 2);
    newTreatImm = newScreen;
    newTreatHpv = newScreen;
    newTreatHyst = newScreen;
    vaxdLmtd = zeros(length(s) - 1 , 1);
    vaxdSchool = vaxdLmtd;
    vaxdCU = vaxdLmtd;
    import java.util.LinkedList
    artDistList = artDistListIn;
    artDist = artDistIn;
    artTreatTracker = zeros(length(s) - 1 , disease , viral , gender , age , risk);
    
    %% Main body of simulation
    for i = 2 : length(s) - 1
        year = currYear + s(i) - 1;
        tspan = [s(i) , s(i + 1)]; % evaluate diff eqs over one time interval
        popIn = popVec(i - 1 , :);
        
        if hpvOn
            % Progression/regression from initial HPV infection to
            % precancer stages and cervical cancer. Differential CC
            % detection by CC stage and HIV status/CD4 count.
            [~ , pop , newCC(i , : , : , :) , ccDeath(i , : , : , :)] ...
                = ode4xtra(@(t , pop) ...
                hpvCCdet(t , pop , hpv_hivClear , rImmuneHiv , c3c2Mults , c2c1Mults , muCC , ...
                normalhpvVaxInds , immunehpvVaxInds , infhpvVaxInds , normalhpvNonVaxInds , ...
                immunehpvNonVaxInds , infhpvNonVaxInds , cin3hpvVaxIndsFrom , ccLochpvVaxIndsTo , ...
                ccLochpvVaxIndsFrom , ccReghpvVaxInds , ccDisthpvVaxInds , ...
                cin3hpvNonVaxIndsFrom , ccLochpvNonVaxIndsTo , ccLochpvNonVaxIndsFrom , ...
                ccReghpvNonVaxInds , ccDisthpvNonVaxInds , cin1hpvVaxInds , ...
                cin2hpvVaxInds , cin3hpvVaxInds , cin1hpvNonVaxInds , ...
                cin2hpvNonVaxInds , cin3hpvNonVaxInds , kInf_Cin1 , kCin1_Cin2 , kCin2_Cin3 , ...
                kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , kCin1_Inf , rNormal_Inf , ...
                rImmune , fImm , kRL , kDR , disease , age , hpvVaxStates , ...
                hpvNonVaxStates , hpvTypeGroups) , tspan , popIn);
            popIn = pop(end , :);  % for next module
            if any(pop(end , :) <  0)
                disp('After hpv')
                break
            end
            
            if (year >= hpvScreenStartYear)
                [dPop , newScreen(i , : , : , : , : , : , : , : , :) , ...
                    newTreatImm(i , : , : , : , : , : , : , : , :) , ...
                    newTreatHpv(i , : , : , : , : , : , : , : , :) , ...
                    newTreatHyst(i , : , : , : , : , : , : , : , :)] ...
                    = hpvScreen(popIn , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , risk , ...
                    screenYrs , screenAlgs , year , stepsPerYear , screenAgeAll , screenAgeS , ...
                    noVaxNoScreen , noVaxToScreen , vaxNoScreen , vaxToScreen , noVaxToScreenTreatImm , ...
                    vaxToScreenTreatImm , noVaxToScreenTreatHpv , vaxToScreenTreatHpv , ...
                    noVaxToScreenTreatVaxHpv , vaxToScreenTreatVaxHpv , noVaxToScreenTreatNonVaxHpv , ...
                    vaxToScreenTreatNonVaxHpv , noVaxToScreenHyst , vaxToScreenHyst , numScreenAge , ageMultsComb , sumall);
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
        [~ , pop , newHpvVax(i , : , : , : , : , :) , newImmHpvVax(i , : , : , : , : , :) , ...
            newHpvNonVax(i , : , : , : , : , :) , newImmHpvNonVax , newHiv(i , : , : , :)] = ...
            ode4xtra(@(t , pop) mixInfect(t , pop , ...
            stepsPerYear , year , disease , intervens , gender , ...
            age , risk , fivYrAgeGrpsOn , hpvTypeGroups , ageSexDebut , gar , epsA_vec , epsR_vec , yr , ...
            partnersM , partnersF , maleActs , femaleActs , ...
            perPartnerHpv_vax , perPartnerHpv_nonV , vaxInds , nonVInds , ...
            lambdaMultImm , lambdaMultVax , artHpvMult , hpv_hivMult , ...
            hpvVaxSus , hpvVaxImm , hpvVaxInf , hpvNonVaxSus , hpvNonVaxImm , hpvNonVaxInf , ...
            circProtect , condProtect , condUse , betaHIVF2M , betaHIVM2F , ...
            hivSus , toHiv , mCurr , fCurr , mCurrArt , fCurrArt , sumall) , tspan , popIn);
        popIn = pop(end , :);
        if any(pop(end , :) < 0)
            disp('After mixInfect')
            break
        end
        
        % HIV module, CD4 Progression, VL progression, ART initiation/dropout,
        % excess HIV mortality
        if hivOn
            [~ , pop , hivDeaths(i , : , :) , artTreat] =...
            ode4xtra(@(t , pop) hiv2a(t , pop , vlAdvancer , artDist , muHIV , ...
            kCD4 ,  maxRateM1 , maxRateF1 , maxRateM2 , maxRateF2 , disease , viral , gender , age , risk , ...
            ageSexDebut , hivInds , stepsPerYear , year , sumall) , tspan , popIn);
            popIn = pop(end , :);    
            artTreatTracker(i , : , : , : , :  ,:) = artTreat;
            artDistList.add(artTreat);
            if artDistList.size() >= stepsPerYear * 2
                artDistList.remove(); % remove CD4 and VL distribution info for people initiating ART more than 2 years ago
            end
            artDist = calcDist(artDistList , disease , viral , gender , age , ...
                risk , sumall); % 2 year average CD4 and VL distribution at time of ART initiation. Details where ART dropouts return to.
            if any(pop(end , :) < 0)
                disp('After hiv')
                break
            end
        end
        
        % Birth, aging, risk redistribution module
        [~ , pop , deaths(i , :)] = ode4xtra(@(t , pop) ...
            bornAgeDieRisk(t , pop , year , ...
            gender , age , fivYrAgeGrpsOn , fertMat , fertMat2 , hivFertPosBirth ,...
        hivFertNegBirth , hivFertPosBirth2 , hivFertNegBirth2 , deathMat , circMat , circMat2 , ...
        MTCTRate , circStartYear , ageInd , riskInd , riskDist , ...
        stepsPerYear , currYear , agesComb , noVaxScreen , noVaxXscreen , ...
        vaxScreen , vaxXscreen , hpvScreenStartYear , sumall) , tspan , popIn);
        popIn = pop(end , :);
        if any(pop(end , :) < 0)
            disp('After bornAgeDieRisk')
            break
        end
        
        if (year >= vaxStartYear)
            % If within first vaxLimitYrs-many vaccine-limited years
            if vaxLimit && ((year - currYear) <= vaxLimitYrs)
                % HPV vaccination module- vaccine limited years
                [dPop , vaxdLmtd(i , :) , vaxRemain] = hpvVaxLmtd(popIn , year , vaxLimitPerYr , ...
                    disease , viral , risk , hpvVaxStates , hpvNonVaxStates , endpoints , ...
                    intervens , vaxCoverL , vaxRemain , vaxGL);
                pop(end , :) = pop(end , :) + dPop;
                popIn = pop(end , :);
                if any(pop(end , :) < 0)
                    disp('After hpvVaxLmtd')
                    break
                end
            
            % If vaccines are not limited
            else
                % HPV vaccination module- school-based vaccination regimen
                [dPop , vaxdSchool(i , :)] = hpvVaxSchool(popIn , disease , viral , risk , ...
                    hpvVaxStates , hpvNonVaxStates , endpoints , intervens , vaxG , vaxAge , ...
                    vaxRate);
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
                        vaxCoverCU , vaxGCU , vaxDiseaseIndsCU);
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
    popLast = popVec(end , :);
    popVec = sparse(popVec); % compress population vectors
    
    filename = ['vaxSimResult' , num2str(simNum)];
    if waning
        filename = ['vaxWaneSimResult' , num2str(simNum)];
    end
    
    parsave(filename , fivYrAgeGrpsOn , tVec ,  popVec , newHiv ,...
        newHpvVax , newImmHpvVax , newHpvNonVax , newImmHpvNonVax , ...
        hivDeaths , deaths , ccDeath , ...
        newCC , artDist , artDistList , artTreatTracker , ... % vaxdLmtd , vaxdSchool , vaxdCU , newScreen , ...
        currYear , lastYear , vaxRate , vaxEff , popLast , pathModifier); % newTreatImm , newTreatHpv , newTreatHyst , ...
end
disp('Done')

%profile viewer

%%
%vaxCEA(pathModifier)

exit(0)
