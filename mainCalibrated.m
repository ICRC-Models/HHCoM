% Main
% Runs simulation over the time period and time step specified by the user.

close all; clear all; clc
% profile clear;

%% Save pre-loaded parameters and pre-calculated indices and matrices
%loadUp2;

%% Load saved parameters
disp('Initializing. Standby...')

paramDir = [pwd , '/Params/'];

% Load saved parameters
load([paramDir, 'generalParams'], 'stepsPerYear' , 'timeStep' , ...
    'disease' , 'viral' , 'hpvVaxStates' , 'hpvNonVaxStates' , 'endpoints' , ...
    'intervens' , 'gender' , 'age' , 'risk' , 'dim' , 'k' , 'toInd' , ...
    'sumall');
load([paramDir, 'demoBehavParams'], 'mInit' , 'fInit' , 'partnersM' , 'partnersF' , ...
    'maleActs' , 'femaleActs' , 'riskDist' , 'fertility' , 'fertility2' , 'mue');
load([paramDir, 'hivParams'], 'betaHIVM2F' , 'betaHIVF2M' , 'muHIV' , 'kVl');
load([paramDir, 'hpvParams'], 'perPartnerHpv' , 'perPartnerHpv_nonV' , 'fImm' , ...
    'kCin1_Inf' , 'kCin2_Cin1' , 'kCin3_Cin2' , 'kCC_Cin3' , ...
    'rNormal_Inf' , 'kInf_Cin1' , 'kCin1_Cin2' , 'kCin2_Cin3' , 'lambdaMultImm');
load([paramDir, 'intervenParams'], 'circ' , 'condUse' , ...
    'maxRateM1' , 'maxRateM2' , 'maxRateF1' , 'maxRateF2' , ...
    'hivStartYear' , 'circStartYear' , 'vaxStartYear' , ...
    'baseline');
load([paramDir , 'calibData'], 'cinPos2014_obs' , 'cinNeg2014_obs' , ...
    'hpv_hiv_obs' , 'hpv_hivNeg_obs' , 'hpv_hivM2008_obs' , 'hpv_hivMNeg2008_obs' , ...
    'hivPrevM_obs' , 'hivPrevF_obs');

%% Load saved indices
load([paramDir,'mixInfectIndices'], 'mCurr' , 'fCurr' , 'mCurrArt' , 'fCurrArt' , ...
    'gar' , 'hivSus' , 'hpvVaxSus' , 'hpvVaxImm' , 'hpvVaxSusVaxd' , 'hpvVaxImmVaxd' , ...
    'hpvNonVaxSus' , 'hpvNonVaxImm' , 'hpvNonVaxSusVaxd' , 'hpvNonVaxImmVaxd' , ...
    'toHiv' , 'hpvVaxInfs' , 'hpvNonVaxInfs')
load([paramDir,'hivIndices'], 'hivInds')
load([paramDir,'hpvIndices'], 'ccLocInds' , 'ccRegInds' , 'ccDistInds' , ...
    'cin1hpvVaxInds' , 'cin2hpvVaxInds' , 'cin3hpvVaxInds' , 'cin1hpvNonVaxInds' , ...
    'cin2hpvNonVaxInds' , 'cin3hpvNonVaxInds' , 'normalhpvVaxInds' , 'immunehpvVaxInds' , ...
    'infhpvVaxInds' , 'normalhpvNonVaxInds' , 'immunehpvNonVaxInds' , 'infhpvNonVaxInds'))
load([paramDir,'ageRiskInds'], 'ageInd' , 'riskInd')

%% Load saved matrices
load([paramDir,'vlAdvancer'], 'vlAdvancer')
load([paramDir,'fertMat'], 'fertMat')
load([paramDir,'hivFertMats'], 'hivFertPosBirth' , 'hivFertNegBirth')
load([paramDir,'fertMat2'], 'fertMat2')
load([paramDir,'hivFertMats2'], 'hivFertPosBirth2' , 'hivFertNegBirth2')
load([paramDir,'deathMat'], 'deathMat')
load([paramDir,'circMat'], 'circMat')
load([paramDir,'circMat2'] , 'circMat2')

%% Model specs
% choose whether to model hysterectomy
hyst = 0;
% choose whether to model HIV
hivOn = 1;
% choose whether to model HPV
hpvOn = 1;
if hpvOn
    disp('HPV module activated.')
end
if hivOn
    disp('HIV module activated')
end

%% Time
startYear = 1910;
c = fix(clock);
currYear = 2020; % c(1); % get the current year
endYear = currYear;
years = endYear - startYear;

%%  Variables/parameters to set based on your scenario

% DIRECTORY TO SAVE RESULTS
pathModifier = 'toNow_091319_singleAge_baseScreen_noBaseVax_nonVhpv'; % ***SET ME***: name for historical run output file 

% VACCINATION
vaxEff = 0.9; % actually bivalent vaccine, but to avoid adding additional compartments, we use nonavalent vaccine and then reduce coverage
waning = 0;    % turn waning on or off

%Parameters for school-based vaccination regimen  % ***SET ME***: coverage for baseline vaccination of 9-year-old girls
vaxAge = [10];
vaxRate = 0.0; %0.86*(0.7/0.9);    % (9 year-old coverage * bivalent vaccine efficacy adjustment)
vaxG = 2;   % indices of genders to vaccinate (1 or 2 or 1,2)

%% Screening
screenYrs = [2000; 2003; 2016; currYear; 2023; 2030; 2045];
hpvScreenStartYear = screenYrs(1);

screenAlgorithm = 1;
screenAlgs{1} = baseline;
screenAlgs{1}.diseaseInds = [1 : disease];

screenAlgs{1}.screenCover_vec = cell(size(screenYrs , 1) - 1, 1); % save data over time interval in a cell array
for i = 1 : size(screenYrs , 1) - 1          % interpolate values at steps within period
    period = [screenYrs(i) , screenYrs(i + 1)];
    screenAlgs{1}.screenCover_vec{i} = interp1(period , screenAlgs{1}.screenCover(i : i + 1 , 1) , ...
        screenYrs(i) : timeStep : screenYrs(i + 1));
end

% Create screening indices
numScreenAge = length(screenAlgs{1}.screenAge);
agesComb = screenAlgs{1}.screenAge;
screenAgeAll = zeros(disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , intervens , length(screenAlgs{1}.screenAge) , risk);
screenAgeS = zeros(disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , 2 , length(screenAlgs{1}.screenAge) , risk);
noVaxNoScreen = zeros(disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , length(screenAlgs{1}.screenAge) , risk);
noVaxToScreen = noVaxNoScreen;
vaxNoScreen = noVaxNoScreen;
vaxToScreen = noVaxNoScreen;
noVaxToScreenTreatImm = zeros(disease , viral , length(screenAlgs{1}.screenAge) , risk);
vaxToScreenTreatImm = noVaxToScreenTreatImm;
noVaxToScreenTreatHpv = noVaxToScreenTreatImm;
vaxToScreenTreatHpv = noVaxToScreenTreatImm;
noVaxToScreenHyst = noVaxToScreenTreatImm;
vaxToScreenHyst = noVaxToScreenTreatImm;
noVaxToScreenTreatVaxHpv = zeros(disease , viral , hpvNonVaxStates , length(screenAlgs{1}.screenAge) , risk);
vaxToScreenTreatVaxHpv = noVaxToScreenTreatVaxHpv;
noVaxToScreenTreatNonVaxHpv = zeros(disease , viral , hpvVaxStates , length(screenAlgs{1}.screenAge) , risk);
vaxToScreenTreatNonVaxHpv = noVaxToScreenTreatNonVaxHpv;
noVaxScreen = zeros(disease*viral*hpvVaxStates*hpvNonVaxStates*endpoints*risk , length(screenAlgs{1}.screenAge));
noVaxXscreen = noVaxScreen;
vaxScreen = noVaxScreen;
vaxXscreen = noVaxScreen;

for aS = 1 : length(screenAlgs{1}.screenAge)
    a = screenAlgs{1}.screenAge(aS);
    
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

%% Initial Population 
MpopStruc = riskDistM;
FpopStruc = riskDistF;

mPop = zeros(age , risk); % distribute initial population size by gender, age risk
fPop = mPop;

if startYear >= 1980;
    for i = 1 : age
        mPop(i , :) = MpopStruc(i, :).* mInit(i) ./ 1.25;
        fPop(i , :) = FpopStruc(i, :).* fInit(i) ./ 1.25;
    end
else
    for i = 1 : age
        mPop(i , :) = MpopStruc(i, :).* mInit(i) ./ (12*1.12);
        fPop(i , :) = FpopStruc(i, :).* fInit(i) ./ (12*1.12);
    end
end

initPop = zeros(dim);
initPop(1 , 1 , 1 , 1 , 1 , 1 , : , :) = mPop; % HIV-, acute infection, HPV Susceptible, no precancer, __, male
initPop(1 , 1 , 1 , 1 , 1 , 2 , : , :) = fPop; % HIV-, acute infection, HPV Susceptible, no precancer, __, female
initPop_0 = initPop;

% Assumes HPV starts before HIV
infected = initPop_0(1 , 1 , 1 , 1 , 1 , : , 16 : 45 , :) * (0.2 * 0.9975); % initial HPV prevalence among age groups 4 - 9 (sexually active) (HIV-)
initPop(1 , 1 , 1 , 1 , 1 , : , 16 : 45 , :) = ...
    initPop_0(1 , 1 , 1 , 1 , 1 , : , 16 : 45 , :) - infected; % moved from HPV

% Omni-HPV type (transition rates weighted by estimated prevalence in population)
initPop(1 , 1 , 2 , 1 , 1 , : , 16 : 45 , :) = 0.5 .* infected; % half moved to vaccine-type HPV+
initPop(1 , 1 , 3 , 1 , 1 , : , 16 : 45 , :) = 0.5 .* infected; % moved to non-vaccine-type HPV+
assert(~any(initPop(:) < 0) , 'Some compartments negative after seeding HPV infections.')

%% Simulation
disp('Start up')
% profile on
disp(' ')

% Initialize time vector
s = 1 : timeStep : years + 1 + timeStep; % stepSize and steps calculated in loadUp.m
% Initialize performance tracking vector
runtimes = zeros(size(s , 2) - 2 , 1);
tVec = linspace(startYear , endYear , length(s) - 1);
% Initialize other vectors
popVec = spalloc(length(s) - 1 , prod(dim) , 10 ^ 8);
popIn = reshape(initPop , prod(dim) , 1); % initial population to "seed" model
popVec(1 , :) = popIn;
newHiv = zeros(length(s) - 1 , gender , age , risk);
newHpv = zeros(length(s) - 1 , gender , disease , age , risk);
newImmHpv = newHpv;
newVaxHpv = newHpv;
newCC = zeros(length(s) - 1 , disease , hpvVaxStates , age);
newCin1 = newCC;
newCin2 = newCC;
newCin3 = newCC;
ccDeath = newCC;
ccTreated = zeros(length(s) - 1 , disease , hpvVaxStates , age , 3); % 3 cancer stages: local, regional, distant
newScreen = zeros(length(s) - 1 , disease , viral , hpvVaxStates , hpvNonVaxStates , numScreenAge , risk , 2);
newTreatImm = newScreen;
newTreatHpv = newScreen;
newTreatHyst = newScreen;
vaxdSchool = zeros(length(s) - 1 , 1);
hivDeaths = zeros(length(s) - 1 , gender , age);
deaths = popVec; 
import java.util.LinkedList
artDistList = LinkedList();
artTreatTracker = zeros(length(s) - 1 , disease , viral , gender , age , risk);
artDist = zeros(disease , viral , gender , age , risk); % initial distribution of inidividuals on ART = 0

%% Main body of simulation
disp(['Simulating period from ' num2str(startYear) ' to ' num2str(endYear) ...
    ' with ' num2str(stepsPerYear), ' steps per year.'])
disp(' ')
disp('Simulation running...')
disp(' ')

% progressbar('Simulation Progress')
for i = 2 : length(s) - 1
    tic
    year = startYear + s(i) - 1;
%     currStep = round(s(i) * stepsPerYear);
%    disp(['current step = ' num2str(startYear + s(i) - 1) ' ('...
%        num2str(length(s) - i) ' time steps remaining until year ' ...
%        num2str(endYear) ')'])
    tspan = [s(i) , s(i + 1)]; % evaluate diff eqs over one time interval
    popIn = popVec(i - 1 , :);
    
    % Add HIV index cases at hivStartYear if not present in initial population
    if (hivOn && (year == hivStartYear))
        % Initialize hiv cases in population
        popIn_init = popIn;
        
        % Create indices
        fromNonHivNonHpv = sort(toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1:gender , 16 : 30 , 1:risk))); 
        toHivNonHpv = sort(toInd(allcomb(3 , 2 , 1 , 1 , 1 , 1:gender , 16 : 30 , 1:risk)));
        fromNonHivHpv = sort(toInd(allcomb(1 , 1 , 2:4 , 1:hpvNonVaxStates , 1 , 1:gender , 16 : 30 , 1:risk))); 
        toHivHpv = sort(toInd(allcomb(3 , 2 , 2:4 , 1:hpvNonVaxStates , 1 , 1:gender , 16 : 30 , 1:risk)));

        % Distribute HIV infections (HPV-)        
        popIn(fromNonHivNonHpv) = (1 - 0.002) .* popIn_init(fromNonHivNonHpv);    % reduce non-HIV infected
        popIn(toHivNonHpv) = (0.002) .* popIn_init(fromNonHivNonHpv);    % increase HIV infected ( male/female, age groups 4-6, med-high risk) (% prevalence)

        % Distribute HIV infections (HPV+)
        popIn(fromNonHivHpv) = (1 - 0.001) .* popIn_init(fromNonHivHpv);    % reduce non-HIV infected
        popIn(toHivHpv) = (0.001) .* popIn_init(fromNonHivHpv);    % increase HIV infected ( male/female, age groups 4-6, med-high risk) (% prevalence)
    end

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
            disease , age , hpvVaxStates , intervens , ...
            rImmuneHiv , hyst , hystInds , hystSusInds , OMEGA) , tspan , popIn);
        popIn = pop(end , :); % for next module
        if any(pop(end , :) <  0)
            disp('After hpv')
            break
        end
        
        if (year >= hpvScreenStartYear)
            [dPop , newScreen(i , : , : , : , : , : , : , :) , ...
                newTreatImm(i , : , : , : , : , : , : , :) , ...
                newTreatHpv(i , : , : , : , : , : , : , :) , ...
                newTreatHyst(i , : , : , : , : , : , : , :)] ...
                = hpvScreen(popIn , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , risk , ...
                    screenYrs , screenAlgs , year , stepsPerYear , screenAgeAll , screenAgeS , ...
                    noVaxNoScreen , noVaxToScreen , vaxNoScreen , vaxToScreen , noVaxToScreenTreatImm , ...
                    vaxToScreenTreatImm , noVaxToScreenTreatHpv , vaxToScreenTreatHpv , ...
                    noVaxToScreenTreatVaxHpv , vaxToScreenTreatVaxHpv , noVaxToScreenTreatNonVaxHpv , ...
                    vaxToScreenTreatNonVaxHpv , noVaxToScreenHyst , vaxToScreenHyst , numScreenAge);
            pop(end , :) = pop(end , :) + dPop;
            popIn = pop(end , :); % for next module
            if any(pop(end , :) <  0)
                disp('After hpv screen')
                break
            end
        end
        %                 [~ , pop] = ode4x(@(t , pop) hpvTreat(t , pop , disease , viral , hpvVaxStates , age , ...
        %                     intervens , detCC , hivCC , muCC , ccRInds , ccSusInds , ...
        %                     hystPopInds , screenFreq , screenCover , hpvSens , ccTreat , ...
        %                     cytoSens , cin1Inds , cin2Inds , cin3Inds ,  normalInds , getHystPopInds ,...
        %                     OMEGA , leep , hystOption , year) , tspan , pop(end , :));
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
        hpvVaxd2NonV , hpvVaxd2NonVScreen , hpvImmVaxd2NonV , hpvImmVaxd2NonVScreen, ...
        hivSus , toHiv , mCurr , fCurr , mCurrArt , fCurrArt , betaHIVF2M , ...
        betaHIVM2F , disease , viral , gender , age , risk , hpvNonVaxStates , hpvVaxStates , ...
        hrInds , nonVInds , intervens , startYear , stepsPerYear , year) , tspan , popIn);
    popIn = pop(end , :); % for next module
    if any(pop(end , :) < 0)
        disp('After mixInfect')
        break
    end
    
    % HIV module, CD4 Progression, VL progression, ART initiation/dropout,
    % excess HIV mortality
    if (hivOn && (year >= hivStartYear))
        [~ , pop , hivDeaths(i , : , :) , artTreat] =...
            ode4xtra(@(t , pop) hiv2a(t , pop , vlAdvancer , artDist , muHIV , ...
            kCD4 ,  maxRateM1 , maxRateM2 , maxRateF1 , maxRateF2 , disease , ...
            viral , gender , age , risk , k , hivInds , ...
            stepsPerYear , year) , tspan , popIn);
        popIn = pop(end , :);
        artTreatTracker(i , : , : , : , :  ,:) = artTreat;
        artDistList.add(artTreat);
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
    
    % Birth, aging, risk redistribution module
    [~ , pop , deaths(i , :)] = ode4xtra(@(t , pop) ...
        bornAgeDieRisk(t , pop , year , ...
        gender , age , fertility , fertMat , fertMat2 , hivFertPosBirth ,...
        hivFertNegBirth , hivFertPosBirth2 , hivFertNegBirth2 , deathMat , circMat , circMat2 , ...
        MTCTRate , circStartYear , ageInd , riskInd , riskDist , ...
        stepsPerYear , currYear , agesComb , noVaxScreen , noVaxXscreen , ...
        vaxScreen , vaxXscreen , hpvScreenStartYear) , tspan , popIn);
    popIn = pop(end , :);
    if any(pop(end , :) < 0)
        disp('After bornAgeDieRisk')
        break
    end 
    
    if ((year >= vaxStartYear) && (vaxRate > 0))
        % HPV vaccination module- school-based vaccination regimen
        [dPop , vaxdSchool(i , :)] = hpvVaxSchool(popIn , k , ...
            disease , viral , risk , hpvVaxStates , hpvNonVaxStates , ...
            intervens , vaxG , vaxAge , vaxRate);
        pop(end , :) = pop(end , :) + dPop;
        if any(pop(end , :) < 0)
            disp('After hpvVaxSchool')
            break
        end
    end

    % add results to population vector
    popVec(i , :) = pop(end , :)';
    runtimes(i) = toc;
    progressbar(i/(length(s) - 1))
end
popLast = popVec(end-1 , :);
disp(['Reached year ' num2str(endYear)])
popVec = sparse(popVec); % compress population vectors

if ~ exist([pwd , '/HHCoM_Results/'])
    mkdir HHCoM_Results
end

savdir = [pwd , '/HHCoM_Results/'];
save(fullfile(savdir , pathModifier) , 'tVec' ,  'popVec' , 'newHiv' ,...
    'newImmHpv' , 'newVaxHpv' , 'newHpv' , 'hivDeaths' , 'deaths' , 'ccDeath' ,... % 'vaxdSchool' , 'newScreen' , 'newTreatImm' , 'newTreatHpv' , 'newTreatHyst' , ...
    'newCC' , 'artDist' , 'artDistList' , ... %artTreatTracker
    'startYear' , 'endYear' , 'popLast');
disp(' ')
disp('Simulation complete.')

% profile viewer

% %% Runtimes
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

% %% Show results
% showResults(pathModifier)

