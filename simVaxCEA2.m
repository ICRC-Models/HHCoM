% Accepts population vector from calibrated natural history model as input

close all;clear all;clc
%profile clear;

% function simVaxCEA2(endYear, vaxCover , vaxEff, vaxAge , waning , varargin)
% p = inputParser;
% addRequired(p , 'endYear');
% addRequired(p , 'vaxCover');
% addRequired(p , 'vaxEff');
% addRequired(p , 'vaxAge');
% addRequired(p , 'waning');
% addOptional(p , 'origEffVec' , []);
% addOptional(p , 't_wane' , 0);
% parse(p , endYear, vaxCover , vaxEff, vaxAge , varargin{:})
% if nargin
%     origEffVec = varargin{1};
%     t_wane = varargin{2};
% end

%%
disp('Start up')

% Load parameters
paramDir = [pwd , '\Params\'];
load([paramDir,'popData'])
load([paramDir,'HIVParams'])
load([paramDir,'general'])
load([paramDir,'mixInfectParams'])
load([paramDir,'vlBeta'])
load([paramDir,'hpvData'])
load([paramDir,'cost_weights'])
% Load indices
load([paramDir,'mixInfectIndices'])
load([paramDir,'hivIndices'])
load([paramDir,'hpvIndices'])
load([paramDir,'hpvTreatIndices'])
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

% Load population
popIn = load([pwd , '\HHCoM_Results\toNow_021319']);
currPop = popIn.popLast;
artDist = popIn.artDist;
artDistList = popIn.artDistList;

% Use calibrated parameters
load([paramDir , 'calibratedParams'])

c = fix(clock);
currYear = c(1); % get the current year
modelYr1 = startYear; % for calculating sexual mixing (this could be streamlined)
stepsPerYear = 6;
timeStep = 1 / stepsPerYear;

%%  Variables/parameters to set based on your scenario

% Directory to save results
pathModifier = 'test_02252019';
if ~ exist([pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\'])
    mkdir ([pwd, '\HHCoM_Results\Vaccine' , pathModifier, '\'])
end

lastYear = 2099; %endYear;
fImm(1 : age) = 1; % all infected individuals who clear HPV get natural immunity

%% Screening
hpvScreen = 1;    % turn HPV DNA testing on or off
dnaTestYrs = [2023; 2030; 2045];
dnaTestCover = [0.45; 0.45; 0.45];

dnaTestCover_vec = cell(size(yr , 1) - 1, 1); % save data over time interval in a cell array
for i = 1 : size(dnaTestYrs , 1) - 1          % interpolate dnaTestCover values at steps within period
    period = [dnaTestYrs(i) , dnaTestYrs(i + 1)];
    dnaTestCover_vec{i} = interp1(period , dnaTestCover(i : i + 1 , 1) , ...
        dnaTestYrs(i) : timeStep : dnaTestYrs(i + 1));
end

% Create screening indices
if hpvScreen    % if present, add indices for screening
    fromNonVaxS = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 1 : 4 , 1 , ... 
        2 , [8,10] , 1 : risk)); 
    toNonVaxS = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 6 , ...
        2 , [8,10] , 1 : risk));
    fromVaxS = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 1 : 4 , 2 , ... 
        2 , [8,10] , 1 : risk)); 
    toVaxS = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 4 , ...
        2 , [8,10] , 1 : risk));
else
    fromNonVaxS = [];    % have to declare these even if hpvScreen=0 because parfor is dumb
    toNonVaxS = [];
    fromVaxS = [];
    toVaxS = [];
end

%% Vaccination

%vaxEff = [0.7 , 0.9];
vaxEff = [0.9];    % used for all vaccine regimens present

%Parameters for school-based vaccination regimen
vaxAge = 3;
vaxCover = [0.8 , 0.9];
vaxG = [2];   % indices of genders to vaccinate (1 or 2 or 1,2)

% Parameters for catch-up vaccination regimen
vaxCU = 0;    % turn catch-up vaccination on or off
vaxAgeCU = [4,5];    % ages catch-up vaccinated
vaxCoverCU = [0.6,0.5];    % coverage for catch-up vaccination by ages catch-up vaccinated
vaxGCU = [2];    % indices of genders to catch-up vaccinate (1 or 2 or 1,2)

% Parameters for vaccination during limited-vaccine years
vaxLimit = 0;    % turn vaccine limit on or off
vaxLimitYrs = 5;    % number years for which vaccines are limited
vaxLimitPerYr = 20000;    % total vaccines available per year for all interventions
vaxAgeL = 5;
vaxCoverL = 0.5;
vaxGL = 2;    % index of gender to vaccinate during limited-vaccine years

% Set up testParams vector for multiple school-based vaccination regimens
%   Example:
%   90% efficacy against 70% of CC types, 100% efficacy against 70% of types, 100% efficacy against 90% of types
%   vaxEff = [0.9 * 0.7 , 0.7 , 0.9]; 
testParams = allcomb(vaxCover , vaxEff); % test scenarios consist of all combinations of school-based vaccine coverage and efficacy
testParams = [testParams ; [0 , 0]]; % Append no vaccine school-based scenario to test scenarios
nTests = size(testParams , 1); % counts number of school-based scenarios to test

if vaxCU
    vaxCoverCUmat = ones(nTests,2) .* vaxCoverCU;
    vaxCoverCUmat(end,:) = 0.0;
else
    vaxCoverCUmat = zeros(nTests,2);    % have to declare these even if vaxCU=0 because parfor is dumb
end
if vaxLimit
    vaxCoverLmat = ones(nTests,1) .* vaxCoverL;
    vaxCoverLmat(end,1) = 0.0;
else
    vaxCoverLmat = zeros(nTests,1);    % have to declare these even if vaxLimit=0 because parfor is dumb
end

% Parameters for waning
waning = 0;    % turn waning on or off

lambdaMultVaxMat = zeros(age , nTests - 1);   % age-based vector for modifying lambda based on vaccination status
vaxEffInd = repmat(1 : length(vaxEff) , 1 , (nTests - 1) /length(vaxEff));
for n = 1 : nTests - 1
    % No waning
    lambdaMultVaxMat(3 : age , n) = vaxEff(vaxEffInd(n));
    
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
        lambdaMultVaxMat(round(effPeriod / 5) + vaxAge - 1 : age , n) = ...
            max(0 , linspace(vaxInit , ...
            vaxInit - kWane * (1 + age - (round(wanePeriod / 5) + vaxAge)) ,...
            age - (round(wanePeriod / 5) + vaxAge) + 2)'); % ensures vaccine efficacy is >= 0
    end
end
lambdaMultVaxMat = [lambdaMultVaxMat , zeros(age , 1)]; % append 0 vaccine protection for no vaccine scenario

% Uncomment the lines below to visualize vaccine efficacy by age group (proxy for waning post-vaccination) 
figure(); plot(lambdaMultVaxMat * 100); 
title('Vaccine Waning'); xlabel('Age Group'); ylabel('Vaccine Protection (%)')
disp(['Simulating period from ' num2str(currYear) ' to ' num2str(lastYear) ...
    ' with ' num2str(stepsPerYear), ' steps per year.'])

% Create vaccine indices
fromNonV = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 1 , ... 
    min(vaxG):max(vaxG) , vaxAge , 1 : risk)); 
toV = toInd(allcomb(1 : disease , 1 : viral , 1 , 9 , 1 , ...
    min(vaxG):max(vaxG) , vaxAge , 1 : risk));
if vaxCU    % if present, add indices for catch-up vaccination regimen
    for aV = 1:length(vaxAgeCU)
        a = vaxAgeCU(aV);
        fromNonVCU(:,aV) = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 1 , ...
            min(vaxGCU):max(vaxGCU) , a , 1 : risk)); 
        toVCU(:,aV) = toInd(allcomb(1 : disease , 1 : viral , 1 , 9 , 1 , ...
            min(vaxGCU):max(vaxGCU) , a , 1 : risk));
    end
else 
    fromNonVCU = [];    % have to declare these even if vaxCU=0 because parfor is dumb
    toVCU = [];
end
if vaxLimit    % if present, add indices for vaccination during limited-vaccine years
    fromNonVL = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 1 , ...
        vaxGL , vaxAgeL , 1 : risk)); 
    toVL = toInd(allcomb(1 : disease , 1 : viral , 1 , 9 , 1 , ...
        vaxGL , vaxAgeL , 1 : risk));
else
    fromNonVL = [];    % have to declare these even if vaxLimit=0 because parfor is dumb
    toVL = [];
end

%% Initialize fixed parameters

% Initialize time vectors
years = lastYear - currYear;
s = 1 : timeStep : years + 1;

% Model specs
hivOn = 1;
hpvOn = 1;

dim = [disease , viral , hpvTypes , hpvStates , periods , gender , age , risk];

% Intervention start years
circStartYear = 1990;
vaxStartYear = currYear;    % (not currently used)

% ART
import java.util.LinkedList
artDistList = popIn.artDistList;
maxRateM_vec = [0.4 , 0.4];    % as of 2013. Scales up from this value in hiv2a. [age 4-6, age >6]
maxRateF_vec = [0.5 , 0.5];    % as of 2013. Scales up from this value in hiv2a. [age 4-6, age >6]
maxRateM1 = maxRateM_vec(1);
maxRateM2 = maxRateM_vec(2);
maxRateF1 = maxRateF_vec(1);
maxRateF2 = maxRateF_vec(2);

%% Run simulation

%profile on
for n = 1 : nTests
    simNum = n;
    vaxEff = testParams(n , 2);
    lambdaMultVax = 1 - lambdaMultVaxMat(: , n);
    vaxRate = testParams(n , 1);
    if vaxCU
        vaxCoverCU = vaxCoverCUmat(n,:);
    end
    if vaxLimit
        vaxRemain = vaxLimitPerYr;
        vaxCoverL = vaxCoverLmat(n);
    end
    % Initialize vectors
    popVec = spalloc(years / timeStep , prod(dim) , 10 ^ 8);
    popIn = currPop; % initial population to "seed" model
    newHiv = zeros(length(s) - 1 , gender , age , risk);
    newHpv = zeros(length(s) - 1 , gender , disease , age , risk);
    newImmHpv = newHpv;
    newVaxHpv = newHpv;
    newCC = zeros(length(s) - 1 , disease , hpvTypes , age);
    ccDeath = newCC;
    ccTreated = zeros(length(s) - 1 , disease , hpvTypes , age , 3); % 3 cancer stages: local, regional, distant
    hivDeaths = zeros(length(s) - 1 , gender , age);
    deaths = zeros(size(popVec));
    vaxd = zeros(length(s) - 1 , 1);
    screend = zeros(length(s) - 1 , 1);
    artTreatTracker = zeros(length(s) - 1 , disease , viral , gender , age , risk);
    popVec(1 , :) = popIn;
    tVec = linspace(currYear , lastYear , size(popVec , 1));
    k = cumprod([disease , viral , hpvTypes , hpvStates , periods , gender , age]);
    artDist = zeros(disease , viral , gender , age , risk); % initial distribution of inidividuals on ART = 0
    %%
    for i = 2 : length(s) - 1
        year = currYear + s(i) - 1;
        currStep = round(s(i) * stepsPerYear);
        tspan = [s(i) , s(i + 1)]; % evaluate diff eqs over one time interval
        popIn = popVec(i - 1 , :);
        
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
                disease , viral , age , hpvTypes , ...
                rImmuneHiv , vaccinated , hystOption) , tspan , popIn);
            popIn = pop(end , :);
            if any(pop(end , :) <  0)
                disp('After hpv')
                break
            end
        end
        
        % HIV and HPV mixing and infection module. Protective effects of condom
        % coverage, circumcision, ART, PrEP (not currently used) are accounted for. 
        [~ , pop , newHpv(i , : , : , : , :) , newImmHpv(i , : , : , : , :) , ...
            newVaxHpv(i , : , : , : , :) , newHiv(i , : , : , :)] = ...
            ode4xtra(@(t , pop) mixInfect(t , pop , currStep , ...
            gar , perPartnerHpv , perPartnerHpv_lr , perPartnerHpv_nonV , ...
            lambdaMultImm , lambdaMultVax , artHpvMult , epsA_vec , epsR_vec , yr , modelYr1 , ...
            circProtect , condProtect , condUse , actsPer , partnersM , partnersF , ...
            hpv_hivMult , hpvSus , hpvImm , toHpv_Imm , hpvVaxd , hpvVaxd2 , toHpv , toHpv_ImmVax , ...
            toHpv_ImmVaxNonV , hivSus , toHiv , mCurr , fCurr , mCurrArt , fCurrArt , ...
            betaHIVF2M , betaHIVM2F , disease , viral , gender , age , risk , hpvStates , hpvTypes , ...
            hrInds , lrInds , hrlrInds , periods , startYear , stepsPerYear , year) , tspan , popIn);
        popIn = pop(end , :); % for next mixing and infection module
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
                stepsPerYear , year) , tspan , pop(end , :));
            artTreatTracker(i , : , : , : , :  ,:) = artTreat;
            artDistList.add(artTreat);
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
        
        % Birth, aging, risk redistribution module
        [~ , pop , deaths(i , :)] = ode4xtra(@(t , pop) ...
            bornAgeDieRisk(t , pop , year , currStep ,...
            gender , age , risk , fertility , fertMat , fertMat2 ,...
            hivFertPosBirth , hivFertNegBirth , hivFertPosBirth2 , ...
            hivFertNegBirth2 , deathMat , circMat , circMat2 , ...
            MTCTRate , circStartYear , ageInd , riskInd , riskDist ,...
            startYear , endYear, stepsPerYear , currYear) , tspan , pop(end , :));
        if any(pop(end , :) < 0)
            disp('After bornAgeDieRisk')
            break
        end
        
        % Screen for HPV
        if hpvScreen
            % Get screening level based on year
            dataYr1 = dnaTestYrs(1);
            dataYrLast = dnaTestYrs(size(dnaTestYrs , 1));
            now = currStep / stepsPerYear + currYear;
            baseYrInd = max(find(now >= dnaTestYrs , 1, 'last') , 1);    % get index of first year <= current year
            baseYr = dnaTestYrs(baseYrInd);
            if currStep < (dataYr1 - currYear) * stepsPerYear    % screening before first year
                screenRate = dnaTestCover_vec{1}(1);
            elseif currStep < (dataYrLast - currYear) * stepsPerYear    % screening between 1st and last year
                screenRate = dnaTestCover_vec{baseYrInd}(currStep - (baseYr - currYear) * stepsPerYear + 1);
            else    % screening after last year
                lastInd = size(dnaTestCover_vec , 1);
                screenRate = dnaTestCover_vec{lastInd}(size(dnaTestCover_vec{lastInd} , 2));
            end
            % Apply screening
            fracScreen = (sum(pop(end , toNonVaxS) , 2) + sum(pop(end , toVaxS) , 2)) / ... % find proportion of population that is currently screened
                    (sum(pop(end , fromNonVaxS) , 2) + sum(pop(end , toNonVaxS) , 2) + sum(pop(end , fromVaxS) , 2) + sum(pop(end , toVaxS) , 2));
            if screenRate - fracScreen > 10 ^ -6 % when proportion screened is below target screening level
                screenCover = max(0 , (screenRate - fracScreen) ./ (1 - fracScreen)); % screen enough people in age group to reach target
                screenGroupNonVax = screenCover .* pop(end , fromNonVaxS);
                screenGroupVax = screenCover .* pop(end , fromVaxS);
                dPop = zeros(size(pop(end , :)));
                dPop(fromNonVaxS) = -screenGroupNonVax;
                dPop(fromVaxS) = -screenGroupVax;
                dPop(toNonVaxS) = screenGroupNonVax;
                dPop(toVaxS) = screenGroupVax;
                pop(end , :) = dPop + pop(end , :); 
                screend(i , :) = screend(i , :) + sumall(screenGroupNonVax) + sumall(screenGroupVax); % count number of people screened at current time step
            end
        end
         
    
    
    
        % Vaccinate for HPV
        
        % If within first vaxLimitYrs-many vaccine-limited years
        if vaxLimit && ((year - currYear) <= vaxLimitYrs)
            if rem(year,1) == 0.0    % reset vaxRemain at start of each new year to the number of available vaccines per year
                vaxRemain = vaxLimitPerYr;
            end
            fracVaxd = sum(pop(end , toVL) , 2) / ... % find proportion of population that is currently vaccinated
                (sum(pop(end , fromNonVL) , 2) + sum(pop(end , toVL) , 2));
            if vaxCoverL - fracVaxd > 10 ^ -6 % when proportion vaccinated is below target vaccination level
                vaxCover = max(0 , (vaxCoverL - fracVaxd) ./ (1 - fracVaxd)); % vaccinate enough people in age group to reach target
                vaxdGroup = vaxCover .* pop(end , fromNonVL);
                if vaxRemain >= 1    % when vaccines remain
                    if (vaxRemain >= sumall(vaxdGroup))    % when remaining vaccines >= targeted coverage
                        usedVax = sumall(vaxdGroup);
                    else    % when remaining vaccines < targeted coverage
                        usedVax = min(sumall(vaxdGroup), vaxRemain);    
                        vaxdGroup = (usedVax .* (pop(end , fromNonVL) ./ sum(pop(end , fromNonVL))))...
                            .* (pop(end , fromNonVL) > 0);    % proportionately divide available vaccines across compartments
                    end
                    vaxRemain = vaxRemain - usedVax;
                else 
                    vaxdGroup = pop(end , fromNonVL) .* 0;    % when vaccines are depleted
                end
                dPop = zeros(size(pop(end , :)));
                dPop(fromNonVL) = -vaxdGroup;
                dPop(toVL) = vaxdGroup;
                pop(end , :) = dPop + pop(end , :); 
                vaxd(i , :) = vaxd(i , :) + sumall(vaxdGroup); % count number of people vaccinated at current time step
            end
        
        % If vaccines are not limited
        else
            % Apply school-based vaccination regimen
            fracVaxd = sum(pop(end , toV) , 2) / ... % find proportion of population that is currently vaccinated
                    (sum(pop(end , fromNonV) , 2) + sum(pop(end , toV) , 2));
            if vaxRate - fracVaxd > 10 ^ -6 % when proportion vaccinated is below target vaccination level
                vaxCover = max(0 , (vaxRate - fracVaxd) ./ (1 - fracVaxd)); % vaccinate enough people in age group to reach target
                vaxdGroup = vaxCover .* pop(end , fromNonV);
                dPop = zeros(size(pop(end , :)));
                dPop(fromNonV) = -vaxdGroup;
                dPop(toV) = vaxdGroup;
                pop(end , :) = dPop + pop(end , :); 
                vaxd(i , :) = vaxd(i , :) + sumall(vaxdGroup); % count number of people vaccinated at current time step
            end
            
            % If present, apply catch-up vaccination regimen
            if vaxCU
                for aV = 1:length(vaxAgeCU)    % apply appropriate coverage rate for each age group catch-up vaccinated
                    fracVaxd = sum(pop(end , toVCU(:,aV)) , 2) / ... % find proportion of population that is currently vaccinated
                        (sum(pop(end , fromNonVCU(:,aV)) , 2) + sum(pop(end , toVCU(:,aV)) , 2));
                    if vaxCoverCU(aV) - fracVaxd > 10 ^ -6 % when proportion vaccinated is below target vaccination level
                        vaxCover = max(0 , (vaxCoverCU(aV) - fracVaxd) ./ (1 - fracVaxd)); % vaccinate enough people in age group to reach target
                        vaxdGroup = vaxCover .* pop(end , fromNonVCU(:,aV));
                        dPop = zeros(size(pop(end , :)));
                        dPop(fromNonVCU(:,aV)) = -vaxdGroup;
                        dPop(toVCU(:,aV)) = vaxdGroup;
                        pop(end , :) = dPop + pop(end , :); 
                        vaxd(i , :) = vaxd(i , :) + sumall(vaxdGroup); % count number of people vaccinated at current time step
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
        newCC , artTreatTracker , vaxd , screend , ccTreated , ...
        currYear , lastYear , vaxRate , vaxEff , popLast , pathModifier);
end
disp('Done')

%profile viewer

%%
vaxCEA(pathModifier)
