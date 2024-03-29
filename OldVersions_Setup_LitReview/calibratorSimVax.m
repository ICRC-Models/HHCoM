% Accepts population vector from calibrated natural history model as input

% Use calibrated parameters
paramDir = [pwd , '/Params/'];
load([paramDir , 'calibratedParams'])

% Load general parameters
paramDir = [pwd , '/Params/'];
load([paramDir, 'general'],'disease','viral','hpvTypes','hpvStates','periods',...
    'gender','age','risk','k','toInd','sumall')
dim = [disease , viral , hpvTypes , hpvStates , periods , gender , age , risk];

% Time
c = fix(clock);
currYear = 2020; %c(1); % get the current year
stepsPerYear = 6;
timeStep = 1 / stepsPerYear;

% Adjust parameters different than calibrated
perPartnerHpv = 0.0045;
condUse = 0.5 * 0.5;
epsA = [0.3 ; 0.3 ; 0.3];
epsR = [0.3 ; 0.3 ; 0.3];
muHIV(11,2) = 0.02;
OMEGA = zeros(age , 1); % hysterectomy rate

% Load indices
paramDir = [pwd , '/Params/'];
load([paramDir,'mixInfectIndices'])
% load([paramDir ,'mixInfectIndices2']) % to load hpvImmVaxd2
load([paramDir,'hivIndices'])
load([paramDir,'hpvIndices'])
load([paramDir,'ageRiskInds'])
% Load matrices
paramDir = [pwd , '/Params/'];
load([paramDir,'ager'])
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

% Use newly calibrated parameters
masterSetMatrix = load([paramDir , 'masterSets_calib_22Aug19_11.dat']); % load most recent parameter sample
paramSet = masterSetMatrix(:,4628);
pIdx = [1,2,5,6,7,8,9,10,19,22,25];
[paramsAll] = genParamStruct();
paramsSub = cell(length(pIdx),1);
startIdx = 1;
for s = 1 : length(pIdx)
    paramsSub{s}.length = paramsAll{pIdx(s)}.length;
    paramsSub{s}.inds = (startIdx : (startIdx + paramsSub{s}.length - 1));
    startIdx = startIdx + paramsSub{s}.length;
end

%% Calibration parameters (FEED FROM LHS)
if any(1 == pIdx)
    idx = find(1 == pIdx);
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
end
if any(2 == pIdx)
    idx = find(2 == pIdx);
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
end
% partnersM(1:2 , 1:risk) = ones(2 , risk) .* 0.00001;
% partnersF(1:2 , 1:risk) = ones(2 , risk) .* 0.00001;
% partnersM(3:age , 1:risk) = [paramSet(1:14) , paramSet(15:28) , paramSet(29:42)];
% partnersF(3:age , 1:risk) = [paramSet(43:56) , paramSet(57:70) , paramSet(71:84)];
% partnersM(10:age , 3) = ones(7 , 1);
% partnersF(10:age , 3) = ones(7 , 1);
if any(5 == pIdx);
    idx = find(5 == pIdx);
    condUse = paramSet(paramsSub{idx}.inds(:));
end
% condUse = paramSet(169);
if any(6 == pIdx);
    idx = find(6 == pIdx);
    %epsA = paramSet(paramsSub{idx}.inds(:));
    epsA = ones(3,1).*paramSet(paramsSub{idx}.inds(:));
end
if any(7 == pIdx);
    idx = find(7 == pIdx);
    %epsR = paramSet(paramsSub{idx}.inds(:));
    epsR = ones(3,1).*paramSet(paramsSub{idx}.inds(:));
end
% epsA = paramSet(170:172);
% epsR = paramSet(173:175);
if any(8 == pIdx)
    idx = find(8 == pIdx);
    maleActsmult = paramSet(paramsSub{idx}.inds(:));
    %rowL = paramsSub{idx}.length/3;
    %rl = paramsSub{idx}.inds(1:rowL);
    %rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    %rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    maleActs(1:2 , 1:risk) = zeros(2 , risk);
    %maleActs(3:age , 1:risk) = [paramSet(rl) , paramSet(rm).*paramSet(rl) , paramSet(rh).*paramSet(rm).*paramSet(rl)];
    maleActs(3:age , 1:risk) = maleActs(3:age , 1:risk) .* maleActsmult;
end
if any(9 == pIdx)
    idx = find(9 == pIdx);
    femaleActsmult = paramSet(paramsSub{idx}.inds(:));
    %rowL = paramsSub{idx}.length/3;
    %rl = paramsSub{idx}.inds(1:rowL);
    %rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    %rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    femaleActs(1:2 , 1:risk) = zeros(2 , risk);
    %femaleActs(3:age , 1: risk) = [paramSet(rl) , paramSet(rm).*paramSet(rl) , paramSet(rh).*paramSet(rm).*paramSet(rl)];
    femaleActs(3:age , 1: risk) = femaleActs(3:age , 1: risk) .* femaleActsmult;
end
% maleActs(1:2 , 1:risk) = zeros(2 , risk);
% femaleActs(1:2 , 1:risk) = zeros(2 , risk);
% maleActs(3:age , 1:risk) = [paramSet(176:189) , paramSet(190:203) , paramSet(204:217)];
% femaleActs(3:age , 1:risk) = [paramSet(218:231) , paramSet(232:245) , paramSet(246:259)];
if any(10 == pIdx)
    idx = find(10 == pIdx);
    perPartnerHpv = paramSet(paramsSub{idx}.inds(:));
end
% perPartnerHpv = paramSet(260);
% perPartnerHpv_lr = paramSet(261);
% perPartnerHpv_nonV = paramSet(262);
% hpv_hivMult = paramSet(263:266) .* 2.0;
% rNormal_Inf_orig = rNormal_Inf;
% rNormal_Inf = paramSet(267) .* rNormal_Inf_orig;
if any(15 == pIdx)
    idx = find(15 == pIdx);
    hpv_hivClear(1,1) = paramSet(paramsSub{idx}.inds(1));
    hpv_hivClear(2,1) = hpv_hivClear(1,1)*paramSet(paramsSub{idx}.inds(2));
    hpv_hivClear(3,1) = hpv_hivClear(2,1)*paramSet(paramsSub{idx}.inds(3));
    hpv_hivClear(4,1) = hpv_hivClear(3,1)*paramSet(paramsSub{idx}.inds(4));
end
% hpv_hivClear = paramSet(283:286) .* 0.5;
if any(16 == pIdx)
    idx = find(16 == pIdx);
    c3c2Mults(4,1) = paramSet(paramsSub{idx}.inds(3));
    c3c2Mults(3,1) = c3c2Mults(4,1)*paramSet(paramsSub{idx}.inds(2));
    c3c2Mults(2,1) = c3c2Mults(3,1)*paramSet(paramsSub{idx}.inds(1));
end
% c3c2Mults = paramSet(287:290) .* 2.0;
if any(17 == pIdx)
    idx = find(17 == pIdx);
    c2c1Mults(4,1) = paramSet(paramsSub{idx}.inds(3));
    c2c1Mults(3,1) = c3c2Mults(4,1)*paramSet(paramsSub{idx}.inds(2));
    c2c1Mults(2,1) = c3c2Mults(3,1)*paramSet(paramsSub{idx}.inds(1));
end
% c2c1Mults = paramSet(291:294) .* 2.0;
% kCCDet = paramSet(295:297);
if any(19 == pIdx)
    idx = find(19 == pIdx);
    lambdaMultImmmult = paramSet(paramsSub{idx}.inds(:));
    lambdaMultImm = lambdaMultImm .* lambdaMultImmmult;
end
% lambdaMultImm = paramSet(298:313);
% maxRateM_vec = paramSet(314:315);
% maxRateF_vec = paramSet(316:317);
if any(22 == pIdx)
    idx = find(22 == pIdx);
    artHpvMult = paramSet(paramsSub{idx}.inds(:));
end
% artHpvMult = paramSet(318) .* 2.0;
if any(25 == pIdx)
    idx = find(25 == pIdx);
    kCCcin3mult = paramSet(paramsSub{idx}.inds(:));
    kCC_Cin3 = kCC_Cin3 .* kCCcin3mult;
end
% kCC_Cin3
% kCD4(1,:,:) = [paramSet(319:323) , paramSet(324:328) , paramSet(329:333) , paramSet(334:338)];
% kCD4(2,:,:) = [paramSet(339:343) , paramSet(344:348) , paramSet(349:353) , paramSet(354:358)];
% kVl(1,:,:) = [paramSet(359:363) , paramSet(364:368) , paramSet(369:373) , paramSet(374:378)];
% kVl(2,:,:) = [paramSet(379:383) , paramSet(384:388) , paramSet(389:393) , paramSet(394:398)];

import java.util.LinkedList
maxRateM_vec = [0.4 , 0.4];    % as of 2013. Scales up from this value in hiv2a. [age 4-6, age >6]
maxRateF_vec = [0.55 , 0.55];    % as of 2013. Scales up from this value in hiv2a. [age 4-6, age >6]
maxRateM1 = maxRateM_vec(1);
maxRateM2 = maxRateM_vec(2);
maxRateF1 = maxRateF_vec(1);
maxRateF2 = maxRateF_vec(2);

epsA_vec = cell(size(yr , 1) - 1, 1); % save data over time interval in a cell array
epsR_vec = cell(size(yr , 1) - 1, 1);
for i = 1 : size(yr , 1) - 1          % interpolate epsA/epsR values at steps within period
    period = [yr(i) , yr(i + 1)];
    epsA_vec{i} = interp1(period , epsA(i : i + 1 , 1) , ...
        yr(i) : timeStep : yr(i + 1));
    epsR_vec{i} = interp1(period , epsR(i : i + 1 , 1) , ...
        yr(i) : timeStep : yr(i + 1));
end

%%  Variables/parameters to set based on your scenario

% LOAD POPULATION
popIn = load([pwd , '/HHCoM_Results/toNow_090319calib_22Aug19_11_4628_2020']); % ***SET ME***: name for historical run input file 
currPop = popIn.popLast;
artDist = popIn.artDist;
artDistList = popIn.artDistList;

% DIRECTORY TO SAVE RESULTS
pathModifier = '090919calib_22Aug19_11_4628_2020_SCE5'; % ***SET ME***: name for simulation output file
if ~ exist([pwd , '/HHCoM_Results/Vaccine' , pathModifier, '/'])
    mkdir ([pwd, '/HHCoM_Results/Vaccine' , pathModifier, '/'])
end

% LAST YEAR & IMMMUNITY
lastYear = 2101; % ***SET ME***: end year of simulation run
fImm(1 : age) = 1; % all infected individuals who clear HPV get natural immunity

% SCREENING
screenAlgorithm = 1; % ***SET ME***: screening algorithm to use (1 for baseline, 2 for CISNET, 3 for WHO)
hivPosScreen = 0; % ***SET ME***: 0 applies same screening algorithm for all HIV states; 1 applies baseline screening to HIV- and selected algorithm for HIV+ 
whoScreenAges = [8 , 10]; % , 6]; % ***SET ME***: [8] for 35, [8,10] for 35&45, [6,8,10] for 25&35&45

% VACCINATION
vaxEff = [0.9];    % 9v-vaccine, used for all vaccine regimens present
waning = 0;    % turn waning on or off

% Parameters for baseline vaccination regimen  % ***SET ME***: coverage for baseline vaccination of 9-year-old girls
vaxAgeB = 2;
vaxCoverB = 0.86; %0.86*(0.7/0.9);    % (9 year-olds: vax whole age group vs. 1/5th (*0.20) to get correct coverage at transition to 10-14 age group) * (bivalent vaccine efficacy adjustment)
vaxGB = 2;   % indices of genders to vaccinate (1 or 2 or 1,2)

%Parameters for school-based vaccination regimen  % ***SET ME***: coverage for school-based vaccination of 10-14 year-old girls
vaxAge = [2 , 3];
vaxCover = [0.86]; %[0.8 , 0.9];
vaxG = [2];   % indices of genders to vaccinate (1 or 2 or 1,2)

% Parameters for catch-up vaccination regimen
vaxCU = 1;    % turn catch-up vaccination on or off  % ***SET ME***: 0 for no catch-up vaccination, 1 for catch-up vaccination
hivPosVaxCU = 0; % ***SET ME***: 0 applies catch-up vaccination algorithm for all HIV states; 1 applies catch-up vaccination only to HIV+ 
vaxAgeCU = [4 , 5 , 6]; %[4:age];   % ages catch-up vaccinated % ***SET ME***: ages for catch-up vaccination
vaxCoverCU = [0.8 , 0.8 , 0.8*0.40]; %ones(1,length(vaxAgeCU)).*0.8;    % coverage for catch-up vaccination by ages catch-up vaccinated % ***SET ME***: coverage for catch-up vaccination by age
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
baseline.screenAge = 8;
baseline.testSens = cytoSens;
baseline.colpoRetain = 0.72;
baseline.cinTreatEff = [0.905 , 0.766 , 0.766 , 0.766 , 0.766 , 0.766 , 0.905 , 0.905 , 0.905 , 0.766]; % cryotherapy/LEEP effectiveness by HIV status
baseline.cinTreatRetain = 0.51;
baseline.cinTreatHpvPersist = 0.28; % HPV persistence with LEEP
baseline.ccTreatRetain = 0.40;
% CISNET
cisnet.screenCover = [0.0; 0.18; 0.48; 0.48; 0.48; 0.70; 0.90];
cisnet.screenAge = [8 , 10];
cisnet.testSens = hpvSens;
cisnet.colpoRetain = 0.81*0.85; % (compliance) * (CIN2+/CC correctly identified by same-day colposcopy)
cisnet.cinTreatEff = baseline.cinTreatEff;
cisnet.cinTreatRetain = 1.0;
cisnet.cinTreatHpvPersist = 0.48; % HPV persistence with cryotherapy 
cisnet.ccTreatRetain = 1.0;
% WHO screening algorithm
who.screenCover = [0.0; 0.18; 0.48; 0.48*0.90; 0.48*0.90; 0.70*0.90; 0.90*0.90]; %90% screening compliance beginning in current year
who.screenAge = whoScreenAges;
who.testSens = hpvSensWHO;
who.colpoRetain = 1.0;
who.cinTreatEff = [1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0];
who.cinTreatRetain = 0.90; % treatment compliance
who.cinTreatHpvPersist = 0.0; %100% treatment efficacy 
who.ccTreatRetain = 0.90; % treatment compliance

if (screenAlgorithm == 1)
    % Baseline screening algorithm
    screenAlgs{1} = baseline;
elseif (screenAlgorithm == 2)
    % CISNET
    screenAlgs{1} = cisnet;
elseif (screenAlgorithm ==3)
    % WHO screening algorithm
    screenAlgs{1} = who;
end

if hivPosScreen
    screenAlgs{2} = baseline;
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
screenAgeAll = zeros(disease , viral , hpvTypes , hpvStates , periods , length(screenAlgs{1}.screenAge) , risk);
screenAgeS = zeros(disease , viral , hpvTypes , hpvStates , 2 , length(screenAlgs{1}.screenAge) , risk);
noVaxNoScreen = zeros(disease , viral , hpvTypes , hpvStates , length(screenAlgs{1}.screenAge) , risk);
noVaxToScreen = noVaxNoScreen;
vaxNoScreen = noVaxNoScreen;
vaxToScreen = noVaxNoScreen;
noVaxToScreenTreatImm = zeros(disease , viral , length(screenAlgs{1}.screenAge) , risk);
vaxToScreenTreatImm = noVaxToScreenTreatImm;
noVaxToScreenTreatHpv = noVaxToScreenTreatImm;
vaxToScreenTreatHpv = noVaxToScreenTreatImm;
noVaxToScreenHyst = noVaxToScreenTreatImm;
vaxToScreenHyst = noVaxToScreenTreatImm;
noVaxScreen = zeros(disease*viral*hpvTypes*hpvStates*risk , length(screenAlgs{1}.screenAge));
noVaxXscreen = noVaxScreen;
vaxScreen = noVaxScreen;
vaxXscreen = noVaxScreen;
for aS = 1 : length(screenAlgs{1}.screenAge)
    a = screenAlgs{1}.screenAge(aS);
    
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
    lambdaMultVaxMat(testParams2{n , 1} : age , n) = vaxEff(vaxEffInd(n));
    
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
        lambdaMultVaxMat(round(effPeriod / 5) + testParams2{n , 1} - 1 : age , n) = ...
            max(0 , linspace(vaxInit , ...
            vaxInit - kWane * (1 + age - (round(wanePeriod / 5) + testParams2{n , 1})) ,...
            age - (round(wanePeriod / 5) + testParams2{n , 1}) + 2)'); % ensures vaccine efficacy is >= 0
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
    ccDeath = newCC;
    ccTreated = zeros(length(s) - 1 , disease , hpvTypes , age , 3); % 3 cancer stages: local, regional, distant
    newScreen = zeros(length(s) - 1 , disease , viral , hpvTypes , hpvStates , risk , 2);
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
                [dPop , newScreen(i , : , : , : , : , : , :) , ...
                    newTreatImm(i , : , : , : , : , : , :) , ...
                    newTreatHpv(i , : , : , : , : , : , :) , ...
                    newTreatHyst(i , : , : , : , : , : , :)] ...
                    = hpvScreen(popIn , disease , viral , hpvTypes , hpvStates , risk , ...
                    screenYrs , screenAlgs , year , stepsPerYear , screenAgeAll , screenAgeS , ...
                    noVaxNoScreen , noVaxToScreen , vaxNoScreen , vaxToScreen , ...
                    noVaxToScreenTreatImm , vaxToScreenTreatImm , noVaxToScreenTreatHpv , ...
                    vaxToScreenTreatHpv , noVaxToScreenHyst , vaxToScreenHyst , ...
                    screenAlgorithm);
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
            gar , perPartnerHpv , perPartnerHpv_lr , perPartnerHpv_nonV , maleActs , femaleActs , ...
            lambdaMultImm , lambdaMultVax , artHpvMult , epsA_vec , epsR_vec , yr , ...
            circProtect , condProtect , condUse , actsPer , partnersM , partnersF , ...
            hpv_hivMult , hpvSus , hpvImm , hpvVaxd , hpvVaxdScreen , hpvVaxd2 , ...
            hpvImmVaxd2 , toHpv , toHpv_Imm , toHpv_Vax , toHpv_VaxScreen , ...
            toHpv_VaxNonV , toHpv_VaxNonVScreen , hivSus , toHiv , mCurr , ...
            fCurr , mCurrArt , fCurrArt , betaHIVF2M , betaHIVM2F , disease , ...
            viral , gender , age , risk , hpvStates , hpvTypes , hrInds , lrInds , ...
            hrlrInds , periods , startYear , stepsPerYear , year) , tspan , popIn);
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
            stepsPerYear , currYear , screenAlgs{1}.screenAge , noVaxScreen , noVaxXscreen , ...
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
                [dPop , vaxdLmtd(i , :) , vaxRemain] = hpvVaxLmtd(popIn , k , ...
                    year , vaxLimitPerYr , disease , viral , risk , hpvTypes , ...
                    hpvStates , periods , vaxCoverL , vaxRemain , vaxGL);
                pop(end , :) = pop(end , :) + dPop;
                popIn = pop(end , :);
                if any(pop(end , :) < 0)
                    disp('After hpvVaxLmtd')
                    break
                end
            
            % If vaccines are not limited
            else
                % HPV vaccination module- school-based vaccination regimen
                [dPop , vaxdSchool(i , :)] = hpvVaxSchool(popIn , k , ...
                    disease , viral , risk , hpvTypes , hpvStates , ...
                    periods , vaxG , vaxAge , vaxRate);
                pop(end , :) = pop(end , :) + dPop;
                popIn = pop(end , :);
                if any(pop(end , :) < 0)
                    disp('After hpvVaxSchool')
                    break
                end
                
                % If present, apply catch-up vaccination regimen
                if vaxCU
                    % HPV vaccination module- catch-up vaccination regimen
                    [dPop , vaxdCU(i , :)] = hpvVaxCU(popIn , k , disease , ...
                        viral , risk , hpvTypes , hpvStates , periods , ...
                        vaxAgeCU , vaxCoverCU , vaxGCU , vaxDiseaseIndsCU);
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
        newCC , artTreatTracker , vaxdLmtd , vaxdSchool , vaxdCU , newScreen , newTreatImm , newTreatHpv , newTreatHyst , ...
        ccTreated , currYear , lastYear , vaxRate , vaxEff , popLast , pathModifier);
end
disp('Done')

%profile viewer

%%
%vaxCEA(pathModifier)
