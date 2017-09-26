function negSumLogL = calibrator(initParams)
hpvOn = 1;
hivOn = 1;
% Load up parameters and data
load('popData')
load('general')
load('general');
load('settings');
load('mixInfectIndices')
load('vlAdvancer')
load('fertMat')
load('hivFertMats')
load('deathMat')
load('circMat')
load('vaxer')
load('mixInfectParams');
load('hpvData')
load('popData')
load('HIVParams')
load('hivIndices')
load('hpvIndices')
load('ager')
load('vlBeta')
load('hpvTreatIndices')
load('calibData')
import java.util.LinkedList
%% Initial population

startYear = 1975;
endYear = 2013; % run to 2013 for calibration
mInit = popInit(: , 1);
MsumInit = sum(mInit);

fInit = popInit(: , 2);
FsumInit = sum(fInit);

MpopStruc = riskDistM;
FpopStruc = riskDistF;

mPop = zeros(age , risk);
fPop = mPop;

for i = 1 : age
    mPop(i , :) = MpopStruc(i, :).* mInit(i) / 1.5;
    fPop(i , :) = FpopStruc(i, :).* fInit(i) / 1.5;
end

dim = [disease , viral , hpvTypes , hpvStates , periods , gender , age ,risk];
initPop = zeros(dim);
initPop(1 , 1 , 1 , 1 , 1 , 1 , : , :) = mPop;
initPop(1 , 1 , 1 , 1 , 1 , 2 , : , :) = fPop;

if hpvOn
    % initPop(1 , 1 , 2 : 4 , 1 , 1 , : , 4 : 9 , :) = 2; % initial HPV hr and lr infecteds (test)
<<<<<<< HEAD
    infected = initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) * 0.1; % try 10% intial HPV prevalence among age groups 4 - 9 (sexually active)
    initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) = ...
        initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) - infected;
    for h = 2
        initPop(1 , 1 , h , 1 , 1 , : , 4 : 9 , :) = infected;
=======
    infected = initPop(1 , 1 , 1 , 1 , 1 , : , 6 : 9 , :) * 0.30; % try 30% intial HPV prevalence among age groups 6 - 9 (sexually active)
    initPop(1 , 1 , 1 , 1 , 1 , : , 6 : 9 , :) = ...
        initPop(1 , 1 , 1 , 1 , 1 , : , 6 : 9 , :) - infected;
    infected45 = initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 5 , :) * 0.60; %try 60% initial HPV prevalence among age groups 4 - 5 (more sexually active)
    initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 5 , :) = ...
        initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 5 , :) - infected45;
    for h = 2
        initPop(1 , 1 , h , 1 , 1 , : , 6 : 9 , :) = infected;
        initPop(1 , 1 , h , 1 , 1 , : , 4 : 5 , :) = infected45;
        initPop(1 , 1 , h , 3 , 1 , : , 4 : 13 , :) = ...
            initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 13 , :) .* 0.07;
        initPop(1 , 1 , h , 4 , 1 , : , 4 : 13 , :) = ...
            initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 13 , :) .* 0.03;
        initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 13 , :) = ...
            initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 13 , :) * 0.9;
>>>>>>> 0b8d00309157244bdeabda53dafa8f57b2fd4f0f
    end
    initPop = max(initPop , 0);
end
assert(~any(initPop(:) < 0) , 'Some compartments negative after seeding HPV infections.')

if hivOn
<<<<<<< HEAD
    toInfectM = (sum(mPop(:)) + sum(fPop(:))) * 0.001* 0.5;
    toInfectF = (sum(mPop(:)) + sum(fPop(:))) * 0.001 * 0.5;
    %     initPop(3 : 4 , 2 : 4 , 1 , 1 , 1 , 1 , 5 , 2) = 1; % initial HIV infected males (acute)
    %     initPop(3 : 4 , 2 : 4 , 1 , 1 , 1 , 2 , 4 , 2) = 1; % initial HIV infected females (acute)
    initPop(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 3) = toInfectM / 3; % initial HIV infected male (0.1% prevalence)
    initPop(1 , 1 , 1 , 1 , 1 , 1 , 4 : 6 , 3) = ...
        initPop(1 , 1 , 1 , 1 , 1 , 1 , 4 : 6 , 3) - toInfectM / 3; % moved to HIV infected
    initPop(3 , 2 , 1 , 1 , 1 , 2 , 4 : 7 , 2) = toInfectF / 6; % initial HIV infected female (0.1% prevalence)
    initPop(1 , 1 , 1 , 1 , 1 , 2 , 4 : 7 , 2) = ...
        initPop(1 , 1 , 1 , 1 , 1 , 2 , 4 : 7 , 2) - toInfectF / 6; % moved to HIV infected
=======
%     toInfectM = (sum(mPop(:)) + sum(fPop(:))) * 0.001* 0.5;
%     toInfectF = (sum(mPop(:)) + sum(fPop(:))) * 0.001 * 0.5;
    %     initPop(3 : 4 , 2 : 4 , 1 , 1 , 1 , 1 , 5 , 2) = 1; % initial HIV infected males (acute)
    %     initPop(3 : 4 , 2 : 4 , 1 , 1 , 1 , 2 , 4 , 2) = 1; % initial HIV infected females (acute)
    initPop(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3) = 0.03 / 2 .* ...
        initPop(1 , 1 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3); % initial HIV infected male (4% prevalence)
    initPop(1 , 1 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3) = ...
        initPop(1 , 1 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3) .* (1 - 0.03 / 2); % moved to HIV infected
    initPop(3 , 2 , 1 , 1 , 1 , 2 , 4 : 6 , 2 : 3) = 0.03 / 2 .*...
        initPop(1 , 1 , 1 , 1 , 1 , 2 , 4 : 6 , 2 : 3); % initial HIV infected female (4% prevalence)
    initPop(1 , 1 , 1 , 1 , 1 , 2 , 4 : 6 , 2 : 3) = ...
        initPop(1 , 1 , 1 , 1 , 1 , 2 , 4 : 6 , 2 : 3) .* (1 - 0.03 / 2); % moved to HIV infected
    
    % HPV infected HIV+
    % females
    initPop(3 , 2 , 2 , 1 , 1 , 2 , 4 : 6 , 1 : 3) = 0.5 .* ...
        initPop(3 , 2 , 1 , 1 , 1 , 2 , 4 : 6 , 1 : 3);
    initPop(3 , 2 , 1 , 1 , 1 , 2 , 4 : 6 , 1 : 3) = 0.5 .* ...
        initPop(3 , 2 , 1 , 1 , 1 , 2 , 4 : 6 , 1 : 3);
    
    % males
    initPop(3 , 2 , 2 , 1 , 1 , 1 , 4 : 6 , 1 : 3) = 0.5 .* ...
        initPop(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 1 : 3);
    initPop(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 1 : 3) = 0.5 .* ...
        initPop(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 1 : 3);   
    
>>>>>>> 0b8d00309157244bdeabda53dafa8f57b2fd4f0f
end
assert(~any(initPop(:) < 0) , 'Some compartments negative after seeding HIV infections.')

% Intervention start years
circStartYear = 1990;
vaxStartYear = 2018;

%% Parameters to be calibrated
fImm(1 : age) = 1; % all individuals get immunity after clearing HPV
<<<<<<< HEAD
lambdaMultImm = initParams(1 : age);
kCin2_Cin3 = initParams(age + 1 : 2 * age);
kCin3_Cin2 = initParams(2 * age + 1 : 3 * age);
kCC_Cin3 = initParams(3 * age + 1 : 4 * age);
perActHpv = initParams(4 * age + 1);
lambdaMultVax = 1 - 0.8 * 0.9;
% analProp = [0 , 0; 
%     0.5111 , 0.4266; 
%     0.5111 , 0.4266]; % risk x gender
% analProp = analProp .* 0; % no anal transmission for now
% analTrans = [138 ; 11] ./ 10 ^ 4; % gender x 1
% vagTransM = 8 / 10 ^ 4 * ones(size(analProp , 1) , 1); 
% vagTransF = 4 / 10 ^ 4 * ones(size(analProp , 1) , 1); 
% transM = vagTransM .* (1 - analProp(: , 1)) + analTrans(1) * analProp(: , 1);
% transF = vagTransF .* (1 - analProp(: , 2)) + analTrans(2) * analProp(: , 2);
% betaHIV_F2M = bsxfun(@times , [7 1 5.8 6.9 11.9 0.04;
%     7 1 5.8 6.9 11.9 0.04;
%     7 1 5.8 6.9 11.9 0.04] , transF);
% betaHIV_M2F = bsxfun(@ times , [7 1 5.8 6.9 11.9 0.04;
%     7 1 5.8 6.9 11.9 0.04;
%     7 1 5.8 6.9 11.9 0.04] , transM);
% 
% yr = [1985; 1988 ; 2003];
% step = 1 / stepsPerYear;
% epsA_vec = cell(size(yr , 1) - 1, 1); % save data over time interval in a cell array
% epsR_vec = cell(size(yr , 1) - 1, 1); 
% for i = 1 : size(yr , 1) - 1
%     period = [yr(i) , yr(i + 1)];
%     epsA_vec{i} = interp1(period , epsA(i : i + 1 , 1) , ...
%         yr(i) : step : yr(i + 1));
%     epsR_vec{i} = interp1(period , epsR(i : i + 1 , 1) , ...
%         yr(i) : step : yr(i + 1));
% end
% 
% betaHIVF2M = zeros(age , risk , viral);
% betaHIVM2F = betaHIVF2M;
% for a = 1 : age
%     betaHIVF2M(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_F2M , maleActs(a , :)')); % HIV(-) males
%     betaHIVM2F(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_M2F , femaleActs(a , :)')); % HIV(-) females
% end
%% Simulation variable preparation
OMEGA(1 : 3) = 0;
OMEGA(4 : age) = logspace(-log(1 - 0.05) , - log(1 - 0.4) , age - 3);
hivPositiveArtAll = sort(toInd(allcomb(10 , 6 , 1 : hpvStates , 1 : hpvTypes , ...
    1 : periods , 1 : gender  , 1 : age , 1 : risk)));
disp(' ')
=======
% lambdaMultImm = initParams(1 : age);
kCin2_Cin3 = initParams(age + 1 : 2 * age);
kCin3_Cin2 = initParams(2 * age + 1 : 3 * age);
kCC_Cin3 = initParams(3 * age + 1 : 4 * age);
perPartnerHpv = initParams(4 * age + 1);
%% Simulation variable preparation
fImm(1 : age) = 1; % all infected individuals who clear HPV get natural immunity
lambdaMultImm(1 : 4) = 1 - 0.01;
lambdaMultImm(5 : 10) = 1 - logspace(log10(0.01) , log10(0.2) , 6);
lambdaMultImm(11 : age ) = lambdaMultImm(10);
lambdaMultVax = 1 - (0.9 * 0.8); % not used for calibration but required as a function input
>>>>>>> 0b8d00309157244bdeabda53dafa8f57b2fd4f0f
% Initialize vectors
timeStep = 1 / stepsPerYear;
years = endYear - startYear;
s = 1 : timeStep : years + 1; % stepSize and steps calculated in loadUp.m
<<<<<<< HEAD
artDistMat = zeros(size(prod(dim) , 20)); % initialize artDistMat to track artDist over past 20 time steps
%performance tracking
runtimes = zeros(size(s , 2) - 2 , 1);
artDistList = LinkedList();
popVec = spalloc(years / timeStep , prod(dim) , 10 ^ 8);
popIn = reshape(initPop , prod(dim) , 1); % initial population to "seed" model
newHiv = zeros(length(s) - 1 , gender , age , risk);
newHpv = newHiv;
newImmHpv = newHpv;
newVaxHpv = newHpv;
newCC = zeros(length(s) - 1 , disease , viral , hpvTypes , age);
=======
%performance tracking
popVec = spalloc(years / timeStep , prod(dim) , 10 ^ 8);
popIn = reshape(initPop , prod(dim) , 1); % initial population to "seed" model
newHiv = zeros(length(s) - 1 , gender , age , risk);
newHpv = zeros(length(s) - 1 , gender , disease , age , risk);
newImmHpv = newHpv;
newVaxHpv = newHpv;
newCC = zeros(length(s) - 1 , disease , viral , hpvTypes , age);
ccDeath = newCC;
>>>>>>> 0b8d00309157244bdeabda53dafa8f57b2fd4f0f
hivDeaths = zeros(length(s) - 1 , age);
deaths = popVec;
artTreatTracker = zeros(length(s) - 1 , disease , viral , gender , age , risk);
popVec(1 , :) = popIn;
tVec = linspace(startYear , endYear , size(popVec , 1));
k = cumprod([disease , viral , hpvTypes , hpvStates , periods , gender , age]);
artDist = zeros(disease , viral , gender , age , risk); % initial distribution of inidividuals on ART = 0

%% Simulation
<<<<<<< HEAD
try
    for i = 2 : length(s) - 1
        year = startYear + s(i) - 1;
        currStep = round(s(i) * stepsPerYear);
        tspan = [s(i) , s(i + 1)]; % evaluate diff eqs over one time interval
        popIn = popVec(i - 1 , :);
        if hivOn
            [~ , pop , newHiv(i , : , : , :)] = ...
=======
for i = 2 : length(s) - 1
    tic
    year = startYear + s(i) - 1;
    currStep = round(s(i) * stepsPerYear);
    disp(['current step = ' num2str(startYear + s(i) - 1) ' ('...
        num2str(length(s) - i) ' time steps remaining until year ' ...
        num2str(endYear) ')'])
    tspan = [s(i) , s(i + 1)]; % evaluate diff eqs over one time interval
    popIn = popVec(i - 1 , :);
    if hivOn
        [~ , pop , newHiv(i , : , : , :)] = ...
>>>>>>> 0b8d00309157244bdeabda53dafa8f57b2fd4f0f
                ode4xtra(@(t , pop) mixInfectHIV(t , pop , currStep , ...
                gar , hivSus , toHiv , mCurr , fCurr , ...
                mCurrArt , fCurrArt ,epsA_vec , epsR_vec , yr , modelYr1 , ...
                circProtect , condProtect , condUse , actsPer , partnersM , partnersF , ...
                betaHIVF2M , betaHIVM2F , disease , viral , gender , age , risk , ...
                hpvStates , hpvTypes , k , periods , stepsPerYear , year) , tspan , popIn);
<<<<<<< HEAD
            popIn = pop(end , :); % for next mixing and infection module
            if any(pop(end , :) < 0)
                disp('After mixInfectHIV')
                break
            end
        end
        
        if hpvOn
            [~ , pop , newHpv(i , : , : , :) , newImmHpv(i , : , : , :) , ...
                newVaxHpv(i , : , : , :)] = ...
                    ode4xtra(@(t , pop) mixInfectHPV(t , pop , currStep , ...
                    gar , perActHpv , lambdaMultImm , lambdaMultVax , epsA_vec , epsR_vec , yr , modelYr1 , ...
                    circProtect , condProtect , condUse , actsPer , partnersM , partnersF , ...
                    hpv_hivMult , hpvSus , hpvImm , hpvVaxd , toHpv , toHpv_ImmVax , ...
                    disease , viral , gender , age , risk , hpvStates , hpvTypes , ...
                    hrInds , lrInds , hrlrInds,  k , periods , stepsPerYear , year) , tspan , popIn);
            if any(pop(end , :) < 0)
                disp('After mixInfectHPV')
                break
            end
        end
        
        if hivOn
            [~ , pop , hivDeaths(i , :) , artTreat] =...
                ode4xtra(@(t , pop) hiv(t , pop , vlAdvancer , artDist , muHIV , ...
                kCD4 , disease , viral , gender , age , risk , k , hivInds , ...
                stepsPerYear , year) , tspan , pop(end , :));
            artTreatTracker(i , : , : , : , :  ,:) = artTreat;
            if any(pop(end , :) < 0)
                disp('After hiv')
                break
            end
            %         [~ , artTreat] = ode4x(@(t , artDist) treatDist(t , popCopy(end , :) , year) , tspan , artDist);
            if size(artDistList) >= 20
                artDistList.remove(); % remove earlier artDist matrix more than "20 time steps old"
            else
                artDistList.add(artTreat);
            end
            artDist = calcDist(artDistList);
        end
        
        if hpvOn
            hystOption = 'on';
            [~ , pop , newCC(i , : , : , : , :)] = ode4xtra(@(t , pop) ...
                hpv(t , pop , immuneInds , infInds , cin1Inds , ...
                cin2Inds , cin3Inds , normalInds , ccInds , ccRegInds , ccDistInds , kInf_Cin1 , kInf_Cin2 , ...
                kCin1_Cin2 , kCin1_Cin3 , kCin2_Cin3 , kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , ...
                kCin1_Inf , kCin2_Inf , kCin3_Cin1 , kNormal_Cin1 , kNormal_Cin2 , ...
                rNormal_Inf , hpv_hivClear , c3c2Mults , c2c1Mults , fImm , kRL , kDR , muCC , ...
                disease , viral , age , hpvTypes , hpvStates , hystOption) , tspan , pop(end , :));
            
            %             [~ , pop] = ode4x(@(t , pop) hpvTreat(t , pop , disease , viral , hpvTypes , age , ...
            %                 periods , detCC , hivCC , muCC , ccRInds , ccSusInds , ...
            %                 hystPopInds , screenFreq , screenCover , hpvSens , ccTreat , ...
            %                 cytoSens , cin1Inds , cin2Inds , cin3Inds ,  normalInds , getHystPopInds ,...
            %                 OMEGA , leep , hystOption , year) , tspan , pop(end , :));
        end
        [~ , pop , deaths(i , :)] = ode4xtra(@(t , pop) bornAgeDie(t , pop , ...
            ager , year , currStep , age , fertility , fertMat , hivFertPosBirth ,...
            hivFertNegBirth , deathMat , circMat , vaxer , MTCTRate , circStartYear , ...
            vaxStartYear , startYear , endYear , stepsPerYear) , tspan , pop(end , :));
        if any(pop(end , :) < 0)
            disp('After bornAgeDie')
            break
        end
        % add results to population vector
        popVec(i , :) = pop(end , :);
    end
    
    popVec = sparse(popVec); % compress population vectors
    negSumLogL = likeFun(popVec , cinPos2008_obs , cinNeg2008_obs ,...
    hpv_hiv_2008_obs , hpv_hivNeg_2008_obs , hivPrevM_obs , hivPrevF_obs , ...
    disease , viral , gender , age , risk , hpvTypes , hpvStates , periods , ...
    startYear , stepsPerYear);
catch
    disp('Simulation aborted.')
    error(['Error after ' num2str(startYear + s(i) - 1) ' ('...
        num2str(length(s) - i) ' time steps remaining until year ' ...
        num2str(endYear) ')'])
    save('errorWorkspace')
    disp('All variables saved to errorWorkspace for debugging.')
    return
=======
        popIn = pop(end , :); % for next mixing and infection module
        if any(pop(end , :) < 0)
            disp('After mixInfectHIV')
            break
        end
    end

    if hpvOn
        [~ , pop , newHpv(i , : , : , : , :) , newImmHpv(i , : , : , : , :) , ...
            newVaxHpv(i , : , : , : , :)] = ...
            ode4xtra(@(t , pop) mixInfectHPV(t , pop , currStep , ...
            gar , perPartnerHpv , lambdaMultImm , lambdaMultVax , epsA_vec , epsR_vec , yr , modelYr1 , ...
            circProtect , condProtect , condUse , actsPer , partnersM , partnersF , ...
            hpv_hivMult , hpvSus , hpvImm , hpvVaxd , toHpv , toHpv_ImmVax , ...
            disease , viral , gender , age , risk , hpvStates , hpvTypes , ...
            hrInds , lrInds , hrlrInds,  k , periods , stepsPerYear , year) , tspan , popIn);
        if any(pop(end , :) < 0)
            disp('After mixInfectHPV')
            break
        end
    end

    if hivOn
        [~ , pop , hivDeaths(i , :) , artTreat] =...
            ode4xtra(@(t , pop) hiv(t , pop , vlAdvancer , artDist , muHIV , ...
            kCD4 , disease , viral , gender , age , risk , k , hivInds , ...
            stepsPerYear , year) , tspan , pop(end , :));
        artTreatTracker(i , : , : , : , :  ,:) = artTreat;
        if any(pop(end , :) < 0)
            disp('After hiv')
            break
        end
        %             [~ , artTreat] = ode4x(@(t , artDist) treatDist(t , popCopy(end , :) , year) , tspan , artDist);
        %             if size(artDistList) >= 20
        %                 artDistList.remove(); % remove earlier artDist matrix more than "20 time steps old"
        %             else
        %                 artDistList.add(artTreat);
        %             end
        %             artDist = calcDist(artDistList);
    end

    if hpvOn
        hystOption = 'on';
        [~ , pop , newCC(i , : , : , : , :) , ccDeath(i , : , : , : , :)] ...
                = ode4xtra(@(t , pop) ...
                hpv(t , pop , immuneInds , infInds , cin1Inds , ...
                cin2Inds , cin3Inds , normalInds , ccInds , ccRegInds , ...
                ccDistInds , kInf_Cin1 , kInf_Cin2 , kCin1_Cin2 , kCin1_Cin3 , ...
                kCin2_Cin3 , kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , kCin1_Inf , ...
                kCin2_Inf , kCin3_Cin1 , kNormal_Cin1 , kNormal_Cin2 , ...
                rNormal_Inf , hpv_hivClear , c3c2Mults , c2c1Mults , fImm ,...
                kRL , kDR , muCC , disease , viral , age , hpvTypes , ...
                hpvStates , k_wane , vaccinated , waned , hystOption) , tspan , pop(end , :));

        %                 [~ , pop] = ode4x(@(t , pop) hpvTreat(t , pop , disease , viral , hpvTypes , age , ...
        %                     periods , detCC , hivCC , muCC , ccRInds , ccSusInds , ...
        %                     hystPopInds , screenFreq , screenCover , hpvSens , ccTreat , ...
        %                     cytoSens , cin1Inds , cin2Inds , cin3Inds ,  normalInds , getHystPopInds ,...
        %                     OMEGA , leep , hystOption , year) , tspan , pop(end , :));
    end
    [~ , pop , deaths(i , :)] = ode4xtra(@(t , pop) bornAgeDie(t , pop , ...
        ager , year , currStep , age , fertility , fertMat , hivFertPosBirth ,...
        hivFertNegBirth , deathMat , circMat , vaxerAger , MTCTRate , circStartYear , ...
        vaxStartYear , vaxRate , startYear , endYear , stepsPerYear) , tspan , pop(end , :));
    if any(pop(end , :) < 0)
        disp('After bornAgeDie')
        break
    end
    % add results to population vector
    popVec(i , :) = pop(end , :)';
    runtimes(i) = toc;
    progressbar(i/(length(s) - 1))
>>>>>>> 0b8d00309157244bdeabda53dafa8f57b2fd4f0f
end