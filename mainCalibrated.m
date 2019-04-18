% Main
% Runs simulation over the time period and time step specified by the
% user.
close all; clear all; clc
profile clear;

% [~ , startYear , endYear , stepsPerYr , IntSteps] = Menu();
% disp('Done with input')

%%
disp('Initializing. Standby...')
disp(' ');

% Use calibrated parameters
paramDir = [pwd , '\Params\'];
load([paramDir , 'calibratedParams'])
perPartnerHpv = 0.0045;
% Load parameters
load([paramDir,'general'])
% Load indices
load([paramDir,'mixInfectIndices'])
% load([paramDir ,'mixInfectIndices2']) % to load hpvImmVaxd2
load([paramDir,'hivIndices'])
load([paramDir,'hpvIndices'])
load([paramDir,'ageRiskInds'])
% Load matrices
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
% choose whether to model hysterectomy (not currently implemented)
hyst = 'off';
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

% Directory to save results
pathModifier = 'toNow_041819';

% Time
c = fix(clock);
currYear = c(1); % get the current year
startYear = 1910; %1980
endYear = currYear;
timeStep = 1 / stepsPerYear;

% Intervention start years
hivStartYear = 1980;
circStartYear = 1990;

% ART
import java.util.LinkedList
artDistList = LinkedList();
maxRateM_vec = [0.40 , 0.40];% maxRateM_arr{sim}; % Maximum ART coverage
maxRateF_vec = [0.55 , 0.55];% maxRateF_arr{sim};
maxRateM1 = maxRateM_vec(1);
maxRateM2 = maxRateM_vec(2);
maxRateF1 = maxRateF_vec(1);
maxRateF2 = maxRateF_vec(2);

%% Screening
screenYrs = [2000; 2003; 2016; currYear; 2023; 2030; 2045];
hpvScreenStartYear = screenYrs(1);
cytoSens = 0.58;
hpvSens = 1.0; % 100 percent efficacy for WHO scenario; 0.89;
cinTreatEff = [0.905 , 0.766 , 0.766 , 0.766 , 0.766 , 0.766 , 0.905 , 0.905 , 0.905 , 0.766]; % cryotherapy effectiveness by HIV status

nhScreenParams = 1;
screenCover = [0.0; 0.18; 0.48; 0.48; 0.48; 0.48; 0.48]; %[0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0];
screenAge = [8];
testSens = cytoSens;
% cryoElig = [1.0 , 0.85 , 0.75 , 0.10 , 0.10 , 0.10];
colpoRetain = 0.72;
cinTreatRetain = 0.51;
cinTreatHpvPersist = 0.28;
ccTreatRetain = 0.40;

screenCover_vec = cell(size(screenYrs , 1) - 1, 1); % save data over time interval in a cell array
for i = 1 : size(screenYrs , 1) - 1          % interpolate dnaTestCover values at steps within period
    period = [screenYrs(i) , screenYrs(i + 1)];
    screenCover_vec{i} = interp1(period , screenCover(i : i + 1 , 1) , ...
        screenYrs(i) : timeStep : screenYrs(i + 1));
end

% Create screening indices
screenAgeAll = zeros(disease , viral , hpvTypes , hpvStates , periods , length(screenAge) , risk);
screenAgeS = zeros(disease , viral , hpvTypes , hpvStates , 2 , length(screenAge) , risk);
noVaxNoScreen = zeros(disease , viral , hpvTypes , hpvStates , length(screenAge) , risk);
noVaxToScreen = noVaxNoScreen;
vaxNoScreen = noVaxNoScreen;
vaxToScreen = noVaxNoScreen;
noVaxToScreenTreatImm = zeros(disease , viral , length(screenAge) , risk);
vaxToScreenTreatImm = noVaxToScreenTreatImm;
noVaxToScreenTreatHpv = noVaxToScreenTreatImm;
vaxToScreenTreatHpv = noVaxToScreenTreatImm;
noVaxToScreenHyst = noVaxToScreenTreatImm;
vaxToScreenHyst = noVaxToScreenTreatImm;
noVaxScreen = zeros(disease*viral*hpvTypes*hpvStates*risk , length(screenAge));
noVaxXscreen = noVaxScreen;
vaxScreen = noVaxScreen;
vaxXscreen = noVaxScreen;
for aS = 1 : length(screenAge)
    a = screenAge(aS);
    
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

%% Initial Population 
mInit = popInit(: , 1); % initial male population size by age
fInit = popInit(: , 2); % initial female population size by age

% use male risk distribution for both genders
riskDistF = riskDistM;

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

dim = [disease , viral , hpvTypes , hpvStates , periods , gender , age ,risk];
initPop = zeros(dim);
initPop(1 , 1 , 1 , 1 , 1 , 1 , : , :) = mPop; % HIV-, acute infection, HPV Susceptible, no precancer, __, male
initPop(1 , 1 , 1 , 1 , 1 , 2 , : , :) = fPop; % HIV-, acute infection, HPV Susceptible, no precancer, __, female
initPop_0 = initPop;

if (hivOn && (hivStartYear == startYear))
    initPop(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3) = 0.005 / 2 .* ...
        initPop_0(1 , 1 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3); % initial HIV infected male (age groups 4-6, med-high risk) (% prevalence)
    initPop(1 , 1 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3) = ...
        initPop_0(1 , 1 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3) .* (1 - 0.005 / 2); % moved to HIV infected
    initPop(3 , 2 , 1 , 1 , 1 , 2 , 4 : 6 , 2 : 3) = 0.005 / 2 .*...
        initPop_0(1 , 1 , 1 , 1 , 1 , 2 , 4 : 6 , 2 : 3); % initial HIV infected female (% prevalence)
    initPop(1 , 1 , 1 , 1 , 1 , 2 , 4 : 6 , 2 : 3) = ...
        initPop_0(1 , 1 , 1 , 1 , 1 , 2 , 4 : 6 , 2 : 3) .* (1 - 0.005 / 2); % moved to HIV infected

    if hpvOn
        initPopHiv_0 = initPop;
        % HIV+ not infected by HPV
        % females
        initPop(3 , 2 , 1 , 1 , 1 , 2 , 4 : 6 , 1 : 3) = 0.3 .* ...
            initPopHiv_0(3 , 2 , 1 , 1 , 1 , 2 , 4 : 6 , 1 : 3);
        % males
        initPop(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 1 : 3) = 0.3 .* ...
            initPopHiv_0(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 1 : 3);

        % HIV+ infected by HPV
        % females
        initPop(3 , 2 , 2 , 1 , 1 , 2 , 4 : 6 , 1 : 3) = 0.7 .* ...
            initPopHiv_0(3 , 2 , 1 , 1 , 1 , 2 , 4 : 6 , 1 : 3);
        % males
        initPop(3 , 2 , 2 , 1 , 1 , 1 , 4 : 6 , 1 : 3) = 0.7 .* ...
            initPopHiv_0(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 1 : 3);
    end
end
assert(~any(initPop(:) < 0) , 'Some compartments negative after seeding HIV infections.')

if (hpvOn && hivOn && (hivStartYear == startYear))
    infected = initPop_0(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) * 0.20; % 20% initial HPV prevalence among age groups 4 - 9 (sexually active) (HIV-)
    initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) = ...
        initPop_0(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) - infected; % moved from HPV-

    % Omni-HPV type (transition rates weighted by estimated prevalence in population)
    initPop(1 , 1 , 2 , 1 , 1 , : , 4 : 9 , :) = infected; % moved to HPV+
end
assert(~any(initPop(:) < 0) , 'Some compartments negative after seeding HPV infections.')

if (hpvOn && ~hivOn) || (hpvOn && hivOn && (hivStartYear > startYear))
    infected = initPop_0(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) * (0.2 * 0.9975); % initial HPV prevalence among age groups 4 - 9 (sexually active) (HIV-)
    initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) = ...
        initPop_0(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) - infected; % moved from HPV

    % Omni-HPV type (transition rates weighted by estimated prevalence in population)
    initPop(1 , 1 , 2 , 1 , 1 , : , 4 : 9 , :) = infected; % moved to HPV+
end
assert(~any(initPop(:) < 0) , 'Some compartments negative after seeding HPV infections.')

%% Simulation
disp('Start up')
profile on
disp(' ')

at = @(x , y) sort(prod(dim)*(y-1) + x);
fImm(1 : age) = 1; % all infected individuals who clear HPV get natural immunity

% Initialize time vectors
years = endYear - startYear;
s = 1 : timeStep : years + 1; % stepSize and steps calculated in loadUp.m

% Initialize performance tracking vector
runtimes = zeros(size(s , 2) - 2 , 1);
% Initialize other vectors
popVec = spalloc(years / timeStep , prod(dim) , 10 ^ 8);
popIn = reshape(initPop , prod(dim) , 1); % initial population to "seed" model
newHiv = zeros(length(s) - 1 , gender , age , risk);
newHpv = zeros(length(s) - 1 , gender , disease , age , risk);
newImmHpv = newHpv;
newVaxHpv = newHpv;
newCC = zeros(length(s) - 1 , disease , hpvTypes , age);
newCin1 = newCC;
newCin2 = newCC;
newCin3 = newCC;
ccDeath = newCC;
ccTreated = zeros(length(s) - 1 , disease , hpvTypes , age , 3); % 3 cancer stages: local, regional, distant
newScreen = zeros(length(s) - 1 , hpvTypes , hpvStates);
vaxd = zeros(length(s) - 1 , 1);
hivDeaths = zeros(length(s) - 1 , gender , age);
deaths = popVec; 
artTreatTracker = zeros(length(s) - 1 , disease , viral , gender , age , risk);
popVec(1 , :) = popIn;
tVec = linspace(startYear , endYear , size(popVec , 1));
k = cumprod([disease , viral , hpvTypes , hpvStates , periods , gender , age]);
artDist = zeros(disease , viral , gender , age , risk); % initial distribution of inidividuals on ART = 0

%% Main body of simulation
disp(['Simulating period from ' num2str(startYear) ' to ' num2str(endYear) ...
    ' with ' num2str(stepsPerYear), ' steps per year.'])
disp(' ')
disp('Simulation running...')
disp(' ')

progressbar('Simulation Progress')
for i = 2 : length(s) - 1
    tic
    year = startYear + s(i) - 1;
%     currStep = round(s(i) * stepsPerYear);
    disp(['current step = ' num2str(startYear + s(i) - 1) ' ('...
        num2str(length(s) - i) ' time steps remaining until year ' ...
        num2str(endYear) ')'])
    tspan = [s(i) , s(i + 1)]; % evaluate diff eqs over one time interval
    popIn = popVec(i - 1 , :);
    
    % Add HIV index cases at hivStartYear if not present in initial population
    if (hivOn && (hivStartYear > startYear) && (year == hivStartYear))
        % Initialize hiv cases in population
        popIn_init = popIn;
        
        % Create indices
        fromNonHivNonHpv = sort(toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1:gender , 4:6 , 1:risk))); 
        toHivNonHpv = sort(toInd(allcomb(3 , 2 , 1 , 1 , 1 , 1:gender , 4:6 , 1:risk)));
        fromNonHivHpv = sort(toInd(allcomb(1 , 1 , 2:4 , 1:hpvStates , 1 , 1:gender , 4:6 , 1:risk))); 
        toHivHpv = sort(toInd(allcomb(3 , 2 , 2:4 , 1:hpvStates , 1 , 1:gender , 4:6 , 1:risk)));

        % Distribute HIV infections (HPV-)        
        popIn(fromNonHivNonHpv) = (1 - 0.002) .* popIn_init(fromNonHivNonHpv);    % reduce non-HIV infected
        popIn(toHivNonHpv) = (0.002) .* popIn_init(fromNonHivNonHpv);    % increase HIV infected ( male/female, age groups 4-6, med-high risk) (% prevalence)

        % Distribute HIV infections (HPV+)
        popIn(fromNonHivHpv) = (1 - 0.001) .* popIn_init(fromNonHivHpv);    % reduce non-HIV infected
        popIn(toHivHpv) = (0.001) .* popIn_init(fromNonHivHpv);    % increase HIV infected ( male/female, age groups 4-6, med-high risk) (% prevalence)
    end

    if hpvOn
        hystOption = 'on';
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
            rImmuneHiv , hystOption) , tspan , popIn);
        popIn = pop(end , :);
        if any(pop(end , :) <  0)
            disp('After hpv')
            break
        end
        
        if (year >= hpvScreenStartYear)
            [~ , pop , newScreen(i , : , :)] ...
                = ode4xtra(@(t , pop) ...
                hpvScreen(t , pop , disease , viral , hpvTypes , hpvStates , risk , ...
                screenYrs , screenCover_vec , testSens , screenAge , ...
                year , stepsPerYear , screenAgeAll , screenAgeS , ...
                noVaxNoScreen , noVaxToScreen , vaxNoScreen , vaxToScreen , ...
                noVaxToScreenTreatImm , vaxToScreenTreatImm , noVaxToScreenTreatHpv , ...
                vaxToScreenTreatHpv , noVaxToScreenHyst , vaxToScreenHyst , ...
                colpoRetain , cinTreatRetain , cinTreatEff , cinTreatHpvPersist , ...
                ccTreatRetain , nhScreenParams) , tspan , popIn);
            popIn = pop(end , :);
            if any(pop(end , :) <  0)
                disp('After hpv screen')
                break
            end
        end
        %                 [~ , pop] = ode4x(@(t , pop) hpvTreat(t , pop , disease , viral , hpvTypes , age , ...
        %                     periods , detCC , hivCC , muCC , ccRInds , ccSusInds , ...
        %                     hystPopInds , screenFreq , screenCover , hpvSens , ccTreat , ...
        %                     cytoSens , cin1Inds , cin2Inds , cin3Inds ,  normalInds , getHystPopInds ,...
        %                     OMEGA , leep , hystOption , year) , tspan , pop(end , :));
    end
    
    % HIV and HPV mixing and infection module. Protective effects of condom
    % coverage, circumcision, ART, PrEP (not currently used) are accounted for. 
    [~ , pop , newHpv(i , : , : , : , :) , newImmHpv(i , : , : , : , :) , ...
        newVaxHpv(i , : , : , : , :) , newHiv(i , : , : , :)] = ...
        ode4xtra(@(t , pop) mixInfect(t , pop , ...
        gar , perPartnerHpv , perPartnerHpv_lr , perPartnerHpv_nonV , maleActs , ...
        femaleActs , lambdaMultImm , lambdaMultVax , artHpvMult , epsA_vec , ...
        epsR_vec , yr , circProtect , condProtect , condUse , actsPer , ...
        partnersM , partnersF , hpv_hivMult , hpvSus , hpvImm , hpvVaxd , ...
        hpvVaxdScreen , hpvVaxd2 , hpvImmVaxd2 , toHpv , toHpv_Imm , toHpv_Vax , ...
        toHpv_VaxScreen , toHpv_VaxNonV , toHpv_VaxNonVScreen , hivSus , toHiv , ...
        mCurr , fCurr , mCurrArt , fCurrArt , betaHIVF2M , betaHIVM2F , disease , ...
        viral , gender , age , risk , hpvStates , hpvTypes , hrInds , lrInds , ...
        hrlrInds , periods , startYear , stepsPerYear , year) , tspan , popIn);
    popIn = pop(end , :); % for next mixing and infection module
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
        stepsPerYear , currYear , screenAge , noVaxScreen , noVaxXscreen , ...
        vaxScreen , vaxXscreen , hpvScreenStartYear) , tspan , popIn);
    if any(pop(end , :) < 0)
        disp('After bornAgeDieRisk')
        break
    end 

    % add results to population vector
    popVec(i , :) = pop(end , :)';
    runtimes(i) = toc;
    progressbar(i/(length(s) - 1))
end
popLast = popVec(end , :);
disp(['Reached year ' num2str(endYear)])
popVec = sparse(popVec); % compress population vectors

if ~ exist([pwd , '\HHCoM_Results\'])
    mkdir HHCoM_Results
end

savdir = [pwd , '\HHCoM_Results\'];%'H:\HHCoM_Results';
save(fullfile(savdir , pathModifier) , 'tVec' ,  'popVec' , 'newHiv' ,...
    'newImmHpv' , 'newVaxHpv' , 'newHpv' , 'hivDeaths' , ...
    'deaths' , 'newCC' , 'artTreatTracker' , 'artDist' , 'artDistList' , ... 
    'startYear' , 'endYear' , 'popLast');
disp(' ')
disp('Simulation complete.')

profile viewer

%% Runtimes
figure()
plot(1 : size(runtimes , 1) , runtimes)
xlabel('Step'); ylabel('Time(s)')
title('Runtimes')
%%
avgRuntime = mean(runtimes); % seconds
stdRuntime = std(runtimes); % seconds
disp(['Total runtime: ' , num2str(sum(runtimes) / 3600) , ' hrs' , ' (' , num2str(sum(runtimes) / 60) , ' mins)']);
disp(['Average runtime per step: ' , num2str(avgRuntime / 60) , ' mins (' , num2str(avgRuntime) , ' secs)']);
disp(['Standard deviation: ' , num2str(stdRuntime / 60) , ' mins (' , num2str(stdRuntime) , ' secs)']);
figure()
h = histogram(runtimes);
title('Runtimes')
ylabel('Frequency')
xlabel('Times (s)')

%% Show results
showResults(pathModifier)
