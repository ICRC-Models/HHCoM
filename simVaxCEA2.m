% Accepts population vector from calibrated natural history model as input

close all;clear all;clc
%profile clear;

%% Load calibrated parameters and reset general parameters based on changes to model
disp('Start up')

% Use calibrated parameters
paramDir = [pwd , '/Params/'];
load([paramDir , 'calibratedParams'])

% Load general parameters and reset changed parameters
paramDir = [pwd , '/Params/'];
load([paramDir, 'general'],'disease','viral','hpvTypes','hpvStates','periods',...
    'gender','age','risk','k')
dim = [disease , viral , hpvTypes , hpvStates , periods , gender , age , risk];
at = @(x , y) sort(prod(dim)*(y-1) + x);
toInd = @(x) (x(: , 8) - 1) * k(7) + (x(: , 7) - 1) * k(6) + (x(: , 6) - 1) * k(5) ...
    + (x(: , 5) - 1) * k(4) + (x(: , 4) - 1) * k(3) + (x(: , 3) - 1) * k(2) ...
    + (x(: , 2) - 1) * k(1) + x(: , 1);
sumall = @(x) sum(x(:));

muHIV(11 , 2) = 0.02;

%% Convert 5-year age groups to 1-year age groups
% Replicate rates across single age groups
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

riskDist = zeros(age , risk , 2);
riskDist(: , : , 1) = riskDistM;
riskDist(: , : , 2) = riskDistF;

% Time
c = fix(clock);
currYear = 2020; % c(1); % get the current year
stepsPerYear = 6;
timeStep = 1 / stepsPerYear;

% Reset  additional changed parameters different than calibrated; load indices and matrices
perPartnerHpv = 0.0045;
condUse = 0.20; %0.5 * 0.5;
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
OMEGA = zeros(age , 1); % hysterectomy rate
betaHIVF2M = zeros(age , risk , viral);
betaHIVM2F = betaHIVF2M;
for a = 1 : age % calculate per-partnership probability of HIV transmission
    betaHIVF2M(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_F2M , maleActs(a , :)')); % HIV(-) males
    betaHIVM2F(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_M2F , femaleActs(a , :)')); % HIV(-) females
end
betaHIVM2F = permute(betaHIVM2F , [2 1 3]); % risk, age, vl
betaHIVF2M = permute(betaHIVF2M , [2 1 3]); % risk, age, vl
% % rNormal_Inf = ones(age,1); % for VCLIR analysis
% % hpv_hivClear = ones(4,1);
% % kCIN1_Inf = zeros(age,1);

% To similate a lower HIV prevalence setting, approximating with decreased per-act HIV transmission probability (80% of initial values)
% % analProp = [0 , 0; 0 , 0; 0 ,0]; % [risk x gender]; proportion practicing anal sex (zero)
% % vagTransM = 8 / 10 ^ 4 * ones(size(analProp , 1) , 1) .* 0.80;
% % vagTransF = 4 / 10 ^ 4 * ones(size(analProp , 1) , 1) .* 0.80;
% % transM = vagTransM .* (1 - analProp(: , 1));
% % transF = vagTransF .* (1 - analProp(: , 2));
% % betaHIV_F2M = bsxfun(@times , [7 1 5.8 6.9 11.9 0.04; 7 1 5.8 6.9 11.9 0.04; 7 1 5.8 6.9 11.9 0.04] , transF);
% % betaHIV_M2F = bsxfun(@times , [7 1 5.8 6.9 11.9 0.04; 7 1 5.8 6.9 11.9 0.04; 7 1 5.8 6.9 11.9 0.04] , transM);
% % betaHIVF2M = zeros(age , risk , viral);
% % betaHIVM2F = betaHIVF2M;
% % for a = 1 : age % calculate per-partnership probability of HIV transmission
% %     betaHIVF2M(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_F2M , maleActs(a , :)')); % HIV(-) males
% %     betaHIVM2F(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_M2F , femaleActs(a , :)')); % HIV(-) females
% % end
% % betaHIVM2F = permute(betaHIVM2F , [2 1 3]); % risk, age, vl
% % betaHIVF2M = permute(betaHIVF2M , [2 1 3]); % risk, age, vl

% Load indices
paramDir = [pwd , '/Params/'];
load([paramDir,'mixInfectIndices'])
% load([paramDir ,'mixInfectIndices2']) % to load hpvImmVaxd2
load([paramDir,'hivIndices'])
load([paramDir,'hpvIndices'])
load([paramDir,'ageRiskInds'])
% Load matrices
paramDir = [pwd , '/Params/'];
%load([paramDir,'ager'])
load([paramDir,'vlAdvancer'])
load([paramDir,'fertMat'])
load([paramDir,'hivFertMats'])
load([paramDir,'fertMat2'])
load([paramDir,'hivFertMats2'])
load([paramDir,'deathMat'])
load([paramDir,'circMat'])
load([paramDir,'circMat2'])

% Model specs
hyst = 0;
hivOn = 1;
hpvOn = 1;

% Intervention start years
circStartYear = 1990;
vaxStartYear = 2014;

% ART
import java.util.LinkedList
maxRateM_vec = [0.4 , 0.4];    % as of 2013. Scales up from this value in hiv2a. [age 4-6, age >6]
maxRateF_vec = [0.55 , 0.55];    % as of 2013. Scales up from this value in hiv2a. [age 4-6, age >6]
maxRateM1 = maxRateM_vec(1);
maxRateM2 = maxRateM_vec(2);
maxRateF1 = maxRateF_vec(1);
maxRateF2 = maxRateF_vec(2);

%%  Variables/parameters to set based on your scenario

% LOAD POPULATION
popIn = load([pwd , '/HHCoM_Results/toNow_091319_singleAge_baseScreen_noBaseVax_2020']); % ***SET ME***: name for historical run input file 
currPop = popIn.popLast;
artDist = popIn.artDist;
artDistList = popIn.artDistList;

% DIRECTORY TO SAVE RESULTS
pathModifier = '091319_WHOP1_SCES12'; % ***SET ME***: name for simulation output file
if ~ exist([pwd , '/HHCoM_Results/Vaccine' , pathModifier, '/'])
    mkdir ([pwd, '/HHCoM_Results/Vaccine' , pathModifier, '/'])
end

% LAST YEAR & IMMMUNITY
lastYear = 2121; % ***SET ME***: end year of simulation run
fImm(1 : age) = 1; % all infected individuals who clear HPV get natural immunity

% SCREENING
screenAlgorithm = 1; % ***SET ME***: screening algorithm to use (1 for baseline, 2 for CISNET, 3 for WHOa, 4 for WHOb)
hivPosScreen = 0; % ***SET ME***: 0 applies same screening algorithm (screenAlgorithm) for all HIV states; 1 applies screenAlgorithm to HIV+ and screenAlgorithmNeg to HIV-
screenAlgorithmNeg = 1; % ***SET ME***: If hivPosScreen=1, screening algorithm to use for HIV- persons (1 for baseline, 2 for CISNET, 3 for WHOa, 4 for WHOb) 
whoScreenAges = [36 , 46]; %[26 , 29 , 32 , 35 , 38 , 41 , 44 , 47 , 50]; % ***SET ME***: ages that get screened when using the WHOa algorithm
cisnetScreenAges = [36 , 46]; % ***SET ME***: ages that get screened when using the CISNET algorithm

% VACCINATION
vaxEff = [0.9];    % 9v-vaccine, used for all vaccine regimens present
waning = 0;    % turn waning on or off

% Parameters for baseline vaccination regimen  % ***SET ME***: coverage for baseline vaccination of 9-year-old girls
vaxAgeB = [10];
vaxCoverB = 0.0; %0.86*(0.7/0.9);    % (9 year-old coverage * bivalent vaccine efficacy adjustment)
vaxGB = 2;   % indices of genders to vaccinate (1 or 2 or 1,2)

%Parameters for school-based vaccination regimen  % ***SET ME***: coverage for school-based vaccination of 9-14 year-old girls
vaxAge = [10 : 15];
vaxCover = [0.8 , 0.9];
vaxG = [2];   % indices of genders to vaccinate (1 or 2 or 1,2)

% Parameters for catch-up vaccination regimen
vaxCU = 0;    % turn catch-up vaccination on or off  % ***SET ME***: 0 for no catch-up vaccination, 1 for catch-up vaccination
hivPosVaxCU = 1; % ***SET ME***: 0 applies catch-up vaccination algorithm for all HIV states; 1 applies catch-up vaccination only to HIV+ 
vaxAgeCU = [16 : 46]; %[16 : 27];    % ages catch-up vaccinated % ***SET ME***: ages for catch-up vaccination
vaxCoverCU = ones(1,length(vaxAgeCU)).*0.80; %0.50 % coverage for catch-up vaccination by ages catch-up vaccinated % ***SET ME***: coverage for catch-up vaccination by age
vaxGCU = [2];    % indices of genders to catch-up vaccinate (1 or 2 or 1,2)

% Parameters for vaccination during limited-vaccine years
vaxLimit = 0;    % turn vaccine limit on or off
vaxLimitYrs = 5;    % number years for which vaccines are limited
vaxLimitPerYr = 20000;    % total vaccines available per year for all interventions
vaxAgeL = 5;
vaxCoverL = 0.5;
vaxGL = 2;    % index of gender to vaccinate during limited-vaccine years

%% Screening
screenYrs = [2000; 2003; 2016; currYear; 2023; 2030; 2045];
hpvScreenStartYear = screenYrs(1);
cytoSens = [0.0 , 0.0 , 0.57 , 0.57 , 0.57 , 0.57 , 0.57 , 0.0 , 0.0 , 0.0]; % pap spear sensitivity by HPV state
hpvSens = [0.0 , 0.0 , 0.881 , 0.881 , 0.881 , 0.881 , 0.881 , 0.0 , 0.0 , 0.0]; % careHPV sensitivity by HPV state
hpvSensWHO = [0.0 , 0.0 , 0.90 , 0.94 , 0.94 , 0.94 , 0.94 , 0.0 , 0.0 , 0.0]; % HPV test sensitivity by HPV state

% Baseline screening algorithm
baseline.screenCover = [0.0; 0.18; 0.48; 0.48; 0.48; 0.48; 0.48];
baseline.screenAge = 36;
baseline.testSens = cytoSens;
baseline.colpoRetain = 0.72;
baseline.cinTreatEff = [0.905 , 0.766 , 0.766 , 0.766 , 0.766 , 0.766 , 0.905 , 0.905 , 0.905 , 0.766]; % cryotherapy/LEEP effectiveness by HIV status
baseline.cinTreatRetain = 0.51;
baseline.cinTreatHpvPersist = 0.28; % HPV persistence with LEEP
baseline.ccTreatRetain = 0.40;
% CISNET
cisnet.screenCover = [0.0; 0.18; 0.48; 0.48; 0.48; 0.70; 0.90];
cisnet.screenAge = cisnetScreenAges; %[36 , 46];
cisnet.testSens = hpvSens;
cisnet.colpoRetain = 0.81*0.85; % (compliance) * (CIN2+/CC correctly identified by same-day colposcopy)
cisnet.cinTreatEff = baseline.cinTreatEff;
cisnet.cinTreatRetain = 1.0;
cisnet.cinTreatHpvPersist = 0.48; % HPV persistence with cryotherapy 
cisnet.ccTreatRetain = 1.0;
% WHO screening algorithm - version a
who.screenCover = [0.0; 0.18; 0.48; 0.48; 0.48; 0.70; 0.90]; % CJB note: removed 90% screening compliance beginning in current year
who.screenAge = whoScreenAges;
who.testSens = hpvSensWHO;
who.colpoRetain = 1.0;
who.cinTreatEff = [1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0];
who.cinTreatRetain = 0.90; % treatment compliance
who.cinTreatHpvPersist = 0.0; %100% treatment efficacy 
who.ccTreatRetain = 0.90; % treatment compliance
% WHO screening algorithm - version b (to apply WHO screening parameters at different ages by HIV status)
whob.screenCover = [0.0; 0.18; 0.48; 0.48; 0.48; 0.70; 0.90]; %CJB note: removed 90% screening compliance beginning in current year
whob.screenAge = [36 , 46];
whob.testSens = hpvSensWHO;
whob.colpoRetain = 1.0;
whob.cinTreatEff = [1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0];
whob.cinTreatRetain = 0.90; % treatment compliance
whob.cinTreatHpvPersist = 0.0; %100% treatment efficacy
whob.ccTreatRetain = 0.90; % treatment compliance

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
    screenAlgs{1}.diseaseInds = [2:6 , 10];
    screenAlgs{2}.diseaseInds = [1 , 7:9];
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
if hivPosScreen
    numScreenAge = numScreenAge + length(screenAlgs{2}.screenAge);
    agesComb = [agesComb , screenAlgs{2}.screenAge];
end
screenAgeAll = zeros(disease , viral , hpvTypes , hpvStates , periods , numScreenAge , risk);
screenAgeS = zeros(disease , viral , hpvTypes , hpvStates , 2 , numScreenAge , risk);
noVaxNoScreen = zeros(disease , viral , hpvTypes , hpvStates , numScreenAge , risk);
noVaxToScreen = noVaxNoScreen;
vaxNoScreen = noVaxNoScreen;
vaxToScreen = noVaxNoScreen;
noVaxToScreenTreatImm = zeros(disease , viral , numScreenAge , risk);
vaxToScreenTreatImm = noVaxToScreenTreatImm;
noVaxToScreenTreatHpv = noVaxToScreenTreatImm;
vaxToScreenTreatHpv = noVaxToScreenTreatImm;
noVaxToScreenHyst = noVaxToScreenTreatImm;
vaxToScreenHyst = noVaxToScreenTreatImm;
noVaxScreen = zeros(disease*viral*hpvTypes*hpvStates*risk , numScreenAge);
noVaxXscreen = noVaxScreen;
vaxScreen = noVaxScreen;
vaxXscreen = noVaxScreen;
for aS = 1 : numScreenAge
    a = agesComb(aS);
    for d = 1 : disease
        for v = 1 : viral
            for h = 1 : hpvTypes
                for s = 1 : hpvStates
                    for r = 1 : risk
                        screenAgeAll(d,v,h,s,:,aS,r) = toInd(allcomb(d , v , h , s , 1 : periods , 2 , a , r)); 
                        screenAgeS(d,v,h,s,:,aS,r) = toInd(allcomb(d , v , h , s , [4,6] , 2 , a , r));

                        noVaxNoScreen(d,v,h,s,aS,r) = sort(toInd(allcomb(d , v , h , s , 1 , 2 , a , r)));
                        noVaxToScreen(d,v,h,s,aS,r) = sort(toInd(allcomb(d , v , h , s , 6 , 2 , a , r)));
                        vaxNoScreen(d,v,h,s,aS,r) = sort(toInd(allcomb(d , v , h , s , 2 , 2 , a , r)));
                        vaxToScreen(d,v,h,s,aS,r) = sort(toInd(allcomb(d , v , h , s , 4 , 2 , a , r)));

                        noVaxToScreenTreatImm(d,v,aS,r) = toInd(allcomb(d , v , 1 , 10 , 6 , 2 , a , r));
                        vaxToScreenTreatImm(d,v,aS,r) = toInd(allcomb(d , v , 1 , 10 , 4 , 2 , a , r));
                        noVaxToScreenTreatHpv(d,v,aS,r) = toInd(allcomb(d , v , 2 , 1 , 6 , 2 , a , r));
                        vaxToScreenTreatHpv(d,v,aS,r) = toInd(allcomb(d , v , 2 , 1 , 4 , 2 , a , r));
                        noVaxToScreenHyst(d,v,aS,r) = toInd(allcomb(d , v , 1 , 8 , 6 , 2 , a , r));
                        vaxToScreenHyst(d,v,aS,r) = toInd(allcomb(d , v , 1 , 8 , 4 , 2 , a , r));
                    end
                end
            end
        end

    end

    % Create indices for removing screening status as people age out of screened age groups
    noVaxScreen(:,aS) = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 6 , ... 
        2 , a+1 , 1 : risk));
    noVaxXscreen(:,aS) = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 , ... 
        2 , a+1 , 1 : risk));
    vaxScreen(:,aS) = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 4 , ... 
        2 , a+1 , 1 : risk));
    vaxXscreen(:,aS) = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 2 , ... 
        2 , a+1 , 1 : risk));
end

%% Vaccination

% Set up differential HIV vaccination for catch-up vaccination regimen
if hivPosVaxCU
    vaxDiseaseIndsCU = [2:6 , 10];
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

lambdaMultVaxMat = zeros(age , nTests); % nTests-1 % age-based vector for modifying lambda based on vaccination status
vaxEffInd = repmat(1 : length(vaxEff) , 1 , (nTests) /length(vaxEff)); % nTests - 1
for n = 1 : nTests %nTests - 1
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
%lambdaMultVaxMat = [lambdaMultVaxMat , zeros(age , 1)]; % append 0 vaccine protection for no vaccine scenario

% Uncomment the lines below to visualize vaccine efficacy by age group (proxy for waning post-vaccination) 
figure(); plot(lambdaMultVaxMat * 100); 
title('Vaccine Waning'); xlabel('Age Group'); ylabel('Vaccine Protection (%)')
disp(['Simulating period from ' num2str(currYear) ' to ' num2str(lastYear) ...
    ' with ' num2str(stepsPerYear), ' steps per year.'])

% Create vaccine indices
% fromNonVSus = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 1 , ... 
%     min(vaxG):max(vaxG) , vaxAge , 1 : risk));
% fromNonVImm = toInd(allcomb(1 : disease , 1 : viral , 2 , 10 , 1 , ... 
%     min(vaxG):max(vaxG) , vaxAge , 1 : risk));
% toV = toInd(allcomb(1 : disease , 1 : viral , 1 , 9 , 1 , ...
%     min(vaxG):max(vaxG) , vaxAge , 1 : risk));
% otherV = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , [1:8,10] , 2:3 , ...
%     min(vaxG):max(vaxG) , vaxAge , 1 : risk));
% allVNonV = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 :hpvStates , 1 : periods , ... 
%     min(vaxG):max(vaxG) , vaxAge , 1 : risk)); 
% if vaxCU    % if present, add indices for catch-up vaccination regimen
%     for aV = 1:length(vaxAgeCU)
%         a = vaxAgeCU(aV);
%         fromNonVCU(:,aV) = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 1 , ...
%             min(vaxGCU):max(vaxGCU) , a , 1 : risk)); 
%         toVCU(:,aV) = toInd(allcomb(1 : disease , 1 : viral , 1 , 9 , 1 , ...
%             min(vaxGCU):max(vaxGCU) , a , 1 : risk));
%     end
% else 
%     fromNonVCU = [];    % have to declare these even if vaxCU=0 because parfor is dumb
%     toVCU = [];
% end
% if vaxLimit    % if present, add indices for vaccination during limited-vaccine years
%     fromNonVL = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 1 , ...
%         vaxGL , vaxAgeL , 1 : risk)); 
%     toVL = toInd(allcomb(1 : disease , 1 : viral , 1 , 9 , 1 , ...
%         vaxGL , vaxAgeL , 1 : risk));
% else
%     fromNonVL = [];    % have to declare these even if vaxLimit=0 because parfor is dumb
%     toVL = [];
% end

%% Run simulation

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
    % Initialize vectors
    years = lastYear - currYear;
    s = 1 : timeStep : years + 1;
    popVec = spalloc(length(s) - 1 , prod(dim) , 10 ^ 8);
    popIn = currPop; % initial population to "seed" model
    newHiv = zeros(length(s) - 1 , gender , age , risk);
    newHpv = zeros(length(s) - 1 , gender , disease , age , risk);
    newImmHpv = newHpv;
    newVaxHpv = newHpv;
    newCC = zeros(length(s) - 1 , disease , hpvTypes , age);
%     newCin1 = newCC;
%     newCin2 = newCC;
%     newCin3 = newCC;
    ccDeath = newCC;
    ccTreated = zeros(length(s) - 1 , disease , hpvTypes , age , 3); % 3 cancer stages: local, regional, distant
    newScreen = zeros(length(s) - 1 , disease , viral , hpvTypes , hpvStates , numScreenAge , risk , 2);
    newTreatImm = newScreen;
    newTreatHpv = newScreen;
    newTreatHyst = newScreen;
    hivDeaths = zeros(length(s) - 1 , gender , age);
    deaths = zeros(size(popVec));
    vaxdLmtd = zeros(length(s) - 1 , 1);
    vaxdSchool = vaxdLmtd;
    vaxdCU = vaxdLmtd;
    artTreatTracker = zeros(length(s) - 1 , disease , viral , gender , age , risk);
    popVec(1 , :) = popIn;
    tVec = linspace(currYear , lastYear-timeStep , length(s)-1);
    k = cumprod([disease , viral , hpvTypes , hpvStates , periods , gender , age]);
    artDist = zeros(disease , viral , gender , age , risk); % initial distribution of inidividuals on ART = 0
    
    %% Main body of simulation
    for i = 2 : length(s) - 1
        year = currYear + s(i) - 1;
        tspan = [s(i) , s(i + 1)]; % evaluate diff eqs over one time interval
        popIn = popVec(i - 1 , :);
        
        if hpvOn
            % Progression/regression from initial HPV infection to
            % precancer stages and cervical cancer. Differential CC
            % detection by CC stage and HIV status/CD4 count.
            [~ , pop , newCC(i , : , : , :) , ccDeath(i , : , : , :) , ...
                ccTreated(i , : , : , : , :)] ...
                = ode4xtra(@(t , pop) ...
                hpvCCdet(t , pop , immuneInds , infInds , cin1Inds , ...
                cin2Inds , cin3Inds , normalInds , ccInds , ccRegInds , ccDistInds , ...
                ccTreatedInds , ccLocDetInds , ccRegDetInds , ccDistDetInds , ...
                kInf_Cin1 , kCin1_Cin2 , kCin2_Cin3 , ...
                kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , kCin1_Inf  ,...
                rNormal_Inf , hpv_hivClear , c3c2Mults , ...
                c2c1Mults , fImm , kRL , kDR , muCC , muCC_det , kCCDet , ...
                disease , age , hpvTypes , ...
                rImmuneHiv , hyst , hystInds , hystSusInds , OMEGA) , tspan , popIn);
            popIn = pop(end , :);  % for next module
            if any(pop(end , :) <  0)
                disp('After hpv')
                break
            end
            
            if (year >= hpvScreenStartYear)
                [dPop , newScreen(i , : , : , : , : , : , : , :) , ...
                    newTreatImm(i , : , : , : , : , : , : , :) , ...
                    newTreatHpv(i , : , : , : , : , : , : , :) , ...
                    newTreatHyst(i , : , : , : , : , : , : , :)] ...
                    = hpvScreen(popIn , disease , viral , hpvTypes , hpvStates , risk , ...
                    screenYrs , screenAlgs , year , stepsPerYear , screenAgeAll , screenAgeS , ...
                    noVaxNoScreen , noVaxToScreen , vaxNoScreen , vaxToScreen , ...
                    noVaxToScreenTreatImm , vaxToScreenTreatImm , noVaxToScreenTreatHpv , ...
                    vaxToScreenTreatHpv , noVaxToScreenHyst , vaxToScreenHyst , ...
                    screenAlgorithm , numScreenAge);
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
        [~ , pop , newHpv(i , : , : , : , :) , newImmHpv(i , : , : , : , :) , ...
            newVaxHpv(i , : , : , : , :) , newHiv(i , : , : , :)] = ...
            ode4xtra(@(t , pop) mixInfect(t , pop , ...
            gar , perPartnerHpv , perPartnerHpv_nonV , maleActs , ...
            femaleActs , lambdaMultImm , lambdaMultVax , artHpvMult , epsA_vec , ...
            epsR_vec , yr , circProtect , condProtect , condUse , actsPer , ...
            partnersM , partnersF , hpv_hivMult , hpvSus , hpvImm , hpvVaxd , ...
            hpvVaxdScreen , hpvVaxd2 , hpvVaxd2Screen , hpvImmVaxd2 , hpvImmVaxd2Screen , ...
            hpvVaxd2NonV , hpvVaxd2NonVScreen , hpvImmVaxd2NonV , hpvImmVaxd2NonVScreen , ...
            hivSus , toHiv , mCurr , fCurr , mCurrArt , fCurrArt , betaHIVF2M , ...
            betaHIVM2F , disease , viral , gender , age , risk , hpvStates , hpvTypes , ...
            hrInds , nonVInds , periods , startYear , stepsPerYear , year) , tspan , popIn);
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
                kCD4 ,  maxRateM1 , maxRateM2 , maxRateF1 , maxRateF2 , disease , ...
                viral , gender , age , risk , k , hivInds , ...
                stepsPerYear , year) , tspan , popIn);
            artTreatTracker(i , : , : , : , :  ,:) = artTreat;
            artDistList.add(artTreat);
            if artDistList.size() >= stepsPerYear * 2
                artDistList.remove(); % remove CD4 and VL distribution info for people initiating ART more than 2 years ago
            end
            artDist = calcDist(artDistList , disease , viral , gender , age , ...
                risk); % 2 year average CD4 and VL distribution at time of ART initiation. Details where ART dropouts return to.
            popIn = pop(end , :);
            if any(pop(end , :) < 0)
                disp('After hiv')
                break
            end
        end
        
        % Birth, aging, risk redistribution module
        [~ , pop , deaths(i , :)] = ode4xtra(@(t , pop) ...
            bornAgeDieRisk(t , pop , year , ...
            gender , age , fertility , fertMat , fertMat2 ,...
            hivFertPosBirth , hivFertNegBirth , hivFertPosBirth2 , ...
            hivFertNegBirth2 , deathMat , circMat , circMat2 , ...
            MTCTRate , circStartYear , ageInd , riskInd , riskDist ,...
            stepsPerYear , currYear , agesComb , noVaxScreen , noVaxXscreen , ...
            vaxScreen , vaxXscreen , hpvScreenStartYear) , tspan , popIn);
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
    
    parsave(filename , tVec ,  popVec , newHiv ,...
        newImmHpv , newVaxHpv , newHpv , deaths , hivDeaths , ccDeath , ...
        newCC , ...%artTreatTracker , vaxdLmtd , vaxdSchool , vaxdCU , newScreen , ccTreated , ...
        currYear , lastYear , vaxRate , vaxEff , popLast , pathModifier); %newTreatImm , newTreatHpv , newTreatHyst , ...
end
disp('Done')

%profile viewer

%%
%vaxCEA(pathModifier)
