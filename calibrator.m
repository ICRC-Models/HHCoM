function negSumLogL = calibrator(initParams)
hpvOn = 1;
hivOn = 1;
% Load up parameters and data
paramDir = [pwd , '\Params\'];
load([paramDir, 'general'])
load([paramDir,'mixInfectIndices'])
load([paramDir,'vlAdvancer'])
load([paramDir,'fertMat'])
load([paramDir,'hivFertMats'])
load([paramDir,'fertMat2'])
load([paramDir,'hivFertMats2'])
load([paramDir,'deathMat'])
load([paramDir,'circMat'])
load([paramDir,'vaxer'])
load([paramDir,'mixInfectParams'])
load([paramDir,'popData'])
load([paramDir,'HIVParams'])
load([paramDir,'hivIndices'])
load([paramDir,'hpvIndices'])
load([paramDir,'ager'])
load([paramDir,'vlBeta'])
load([paramDir,'hpvTreatIndices'])
load([paramDir,'calibParams'])
load([paramDir,'vaxInds'])
load([paramDir,'settings'])
load([paramDir,'hpvData'])
load([paramDir,'calibData'])
load([paramDir,'calibInitParams'])
load([paramDir ,'cost_weights'])
import java.util.LinkedList
vaxerAger = ager;
vaxRate = 0;
startYear = 1975;
endYear = 2015; % run to 2015 for calibration
%% Initial population
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
initPop_0 = initPop;
if hivOn
    %     toInfectM = (sum(mPop(:)) + sum(fPop(:))) * 0.001* 0.5;
    %     toInfectF = (sum(mPop(:)) + sum(fPop(:))) * 0.001 * 0.5;
    %     initPop(3 : 4 , 2 : 4 , 1 , 1 , 1 , 1 , 5 , 2) = 1; % initial HIV infected males (acute)
    %     initPop(3 : 4 , 2 : 4 , 1 , 1 , 1 , 2 , 4 , 2) = 1; % initial HIV infected females (acute)
    initPop(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3) = 0.005 / 2 .* ...
        initPop_0(1 , 1 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3); % initial HIV infected male (% prevalence)
    initPop(1 , 1 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3) = ...
        initPop_0(1 , 1 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3) .* (1 - 0.005 / 2); % moved to HIV infected
    initPop(3 , 2 , 1 , 1 , 1 , 2 , 4 : 6 , 2 : 3) = 0.005 / 2 .*...
        initPop_0(1 , 1 , 1 , 1 , 1 , 2 , 4 : 6 , 2 : 3); % initial HIV infected female (% prevalence)
    initPop(1 , 1 , 1 , 1 , 1 , 2 , 4 : 6 , 2 : 3) = ...
        initPop_0(1 , 1 , 1 , 1 , 1 , 2 , 4 : 6 , 2 : 3) .* (1 - 0.005 / 2); % moved to HIV infected
    
        if hpvOn
            initPopHiv_0 = initPop;
            % HPV infected HIV+
            % females
            initPop(3 , 2 , 1 , 1 , 1 , 2 , 4 : 6 , 1 : 3) = 0.3 .* ...
                initPopHiv_0(3 , 2 , 1 , 1 , 1 , 2 , 4 : 6 , 1 : 3);
    
            % males
            initPop(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 1 : 3) = 0.3 .* ...
                initPopHiv_0(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 1 : 3);
    
            for h = 2 : 4
                % females
                initPop(3 , 2 , h , 1 , 1 , 2 , 4 : 6 , 1 : 3) = 0.7 / 3 .* ...
                    initPopHiv_0(3 , 2 , 1 , 1 , 1 , 2 , 4 : 6 , 1 : 3);
               % males
                initPop(3 , 2 , h , 1 , 1 , 1 , 4 : 6 , 1 : 3) = 0.7 / 3 .* ...
                    initPopHiv_0(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 1 : 3);
            end
        end
end
assert(~any(initPop(:) < 0) , 'Some compartments negative after seeding HIV infections.')

if hpvOn
    infected = initPop_0(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) * 0.20; % 20% intial HPV prevalence among age groups 4 - 9 (sexually active)
    initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) = ...
        initPop_0(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) - infected;

    % HPV 16/18
    initPop(1 , 1 , 2 , 1 , 1 , : , 4 : 9 , :) = 0.7 .* infected;
    
    % 4v and oHR
    for h = 3 : 4
        for s = 2 : 4
            initPop(1 , 1 , h , s , 1 , : , 4 : 9 , :) = 0.3 / 6 * infected; % in CIN1 - CIN 3 states
        end
    end
%     initPop = max(initPop , 0);
end
assert(~any(initPop(:) < 0) , 'Some compartments negative after seeding HPV infections.')

% Intervention start years
circStartYear = 1990;
vaxStartYear = 2017;
%% calibration parameters
% for calibrateModel
% kCin2_Cin3(: , 1) = initParams(1 : age);
% kCin3_Cin2(: , 1) = initParams(age + 1 : 2 * age);
% kCC_Cin3(: , 1) = initParams(2 * age + 1 : 3 * age);
% kCin2_Cin3(: , 2) = initParams(3 * age + 1 : 4 * age);
% kCin3_Cin2(: , 2) = initParams(4 * age + 1 : 5 * age);
% kCC_Cin3(: , 2) = initParams(5 * age + 1 : 6 * age);
% kCin2_Cin3(: , 3) = initParams(6 * age + 1 : 7 * age);
% kCin3_Cin2(: , 3) = initParams(7 * age + 1 : 8 * age);
% kCC_Cin3(: , 3) = initParams(8 * age + 1 : 9 * age);
% rImmuneHiv = initParams(9 * age + 1 : 9 * age + 1 + 3);
% c3c2Mults = initParams(9 * age + 5 : 9 * age + 8);
% c2c1Mults = initParams(9 * age + 9 : 9 * age + 12);
% artHpvMult = initParams(9 * age + 13);
% perPartnerHpv = initParams(9 * age + 14);

% for quickCalibrateModel
for i = 1 : 3
    kCin1_Inf(: , i) = initParams(i) .* kCin1_Inf(: , i);
    rNormal_Inf(: , i) = initParams(3 + i) .* rNormal_Inf(: , i);
    kCC_Cin3(: , i) = initParams(6 + i) .* kCC_Cin3(: , i);
    kCin3_Cin2(: , i) = initParams(49 + i) .* kCin3_Cin2(: , i);
end

rImmuneHiv = initParams(10 : 13);
c3c2Mults = initParams(14 : 17);
c2c1Mults = initParams(18 : 21);
artHpvMult = initParams(22);
perPartnerHpv= initParams(23);
lambdaMultImm = initParams(24 : 39);
hpv_hivClear = initParams(40 : 43);
hpvClearMult = initParams(44 : 47);
perPartnerHpv_lr = initParams(48);
perPartnerHpv_nonV = initParams(49);

for i = 0 : 2
    maleActs(: , i + 1) = maleActs(: , i + 1) .* initParams(53 + i);
    femaleActs(: , i + 1) = femaleActs(: , i + 1) .* initParams(56 + i);
end

for a = 1 : age
    betaHIVF2M(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_F2M , maleActs(a , :)')); % HIV(-) males
    betaHIVM2F(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_M2F , femaleActs(a , :)')); % HIV(-) females
end
betaHIVM2F = permute(betaHIVM2F , [2 1 3]); % risk, age, vl
betaHIVF2M = permute(betaHIVF2M , [2 1 3]); % risk, age, vl
%% Simulation variable preparation
OMEGA(1 : 3) = 0;
OMEGA(4 : age) = logspace(-log(1 - 0.05) , - log(1 - 0.4) , age - 3);
hivPositiveArtAll = sort(toInd(allcomb(10 , 6 , 1 : hpvStates , 1 : hpvTypes , ...
    1 : periods , 1 : gender  , 1 : age , 1 : risk)));
fImm(1 : age) = 1; % all infected individuals who clear HPV get natural immunity
% lambdaMultImm(1 : 4) = 1 - 0.01;
% lambdaMultImm(5 : 10) = 1 - logspace(log10(0.01) , log10(0.1) , 6);
% lambdaMultImm(11 : age) = lambdaMultImm(10);
lambdaMultVax = ones(age , 2);
% Initialize vectors
timeStep = 1 / stepsPerYear;
years = endYear - startYear;
s = 1 : timeStep : years + 1; % stepSize and steps calculated in loadUp.m
artDistMat = zeros(size(prod(dim) , 20)); % initialize artDistMat to track artDist over past 20 time steps
popVec = spalloc(years / timeStep , prod(dim) , 10 ^ 8);
popIn = reshape(initPop , prod(dim) , 1); % initial population to "seed" model
newHiv = zeros(length(s) - 1 , gender , age , risk);
newHpv = zeros(length(s) - 1 , gender , disease , age , risk);
newImmHpv = newHpv;
newVaxHpv = newHpv;
newCC = zeros(length(s) - 1 , disease , hpvTypes , age);
ccDeath = newCC;
ccTreatCost = zeros(length(s) - 1 , disease , hpvTypes , age , 3); % 3 cancer stages: local, regional, distant
vaxd = zeros(length(s) - 1 , 1);
hivDeaths = zeros(length(s) - 1 , gender , age);
deaths = popVec;
artTreatTracker = zeros(length(s) - 1 , disease , viral , gender , age , risk);
popVec(1 , :) = popIn;
k = cumprod([disease , viral , hpvTypes , hpvStates , periods , gender , age]);
artDist = zeros(disease , viral , gender , age , risk); % initial distribution of inidividuals on ART = 0
vaxMat = ager .* 0;

maxRateM_vec = [0.45 , 0.45];% maxRateM_arr{sim};
maxRateF_vec = [0.65 , 0.65];% maxRateF_arr{sim};

maxRateM1 = 1 - exp(-maxRateM_vec(1));
maxRateM2 = 1 - exp(-maxRateM_vec(2));
maxRateF1 = 1 - exp(-maxRateF_vec(1));
maxRateF2 = 1 - exp(-maxRateF_vec(2));
%% Simulation
for i = 2 : length(s) - 1
    tic
    year = startYear + s(i) - 1;
    currStep = round(s(i) * stepsPerYear);
    tspan = [s(i) , s(i + 1)]; % evaluate diff eqs over one time interval
    popIn = popVec(i - 1 , :);
        
    if hpvOn
        hystOption = 'on';
        [~ , pop , newCC(i , : , : , :) , ccDeath(i , : , : , :) , ...
            ccTreatCost(i , : , : , : , :)] ...
            = ode4xtra(@(t , pop) ...
            hpv(t , pop , immuneInds , infInds , cin1Inds , ...
            cin2Inds , cin3Inds , normalInds , ccInds , ccRegInds , ccDistInds , ...
            ccTreatedInds , kInf_Cin1 , kCin1_Cin2 , kCin2_Cin3 , ...
            kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , kCin1_Inf  ,...
            rNormal_Inf , hpv_hivClear , c3c2Mults , hpvClearMult , ...
            c2c1Mults , fImm , kRL , kDR , muCC , kCCDet , ...
            disease , viral , age , hpvTypes , ...
            rImmuneHiv , vaccinated , hystOption) , tspan , popIn);
        popIn = pop(end , :);
        if any(pop(end , :) <  0)
            disp('After hpv')
            break
        end
        
    end
    
    [~ , pop , newHpv(i , : , : , : , :) , newImmHpv(i , : , : , : , :) , ...
        newVaxHpv(i , : , : , : , :) , newHiv(i , : , : , :)] = ...
        ode4xtra(@(t , pop) mixInfect(t , pop , currStep , ...
        gar , perPartnerHpv , perPartnerHpv_lr , perPartnerHpv_nonV , ...
        lambdaMultImm , lambdaMultVax , artHpvMult , epsA_vec , epsR_vec , yr , modelYr1 , ...
        circProtect , condProtect , condUse , actsPer , partnersM , partnersF , ...
        hpv_hivMult , hpvSus , hpvImm , toHpv_Imm , hpvVaxd , hpvVaxd2 , toHpv , toHpv_ImmVax , ...
        hivSus , toHiv , mCurr , fCurr , mCurrArt , fCurrArt , ...
        betaHIVF2M , betaHIVM2F , disease , viral , gender , age , risk , hpvStates , hpvTypes , ...
        hrInds , lrInds , hrlrInds , periods , startYear , stepsPerYear , year) , tspan , popIn);
    popIn = pop(end , :); % for next mixing and infection module
    if any(pop(end , :) < 0)
        disp('After mixInfect')
        break
    end
    
    if hivOn
        [~ , pop , hivDeaths(i , : , :) , artTreat] =...
            ode4xtra(@(t , pop) hiv2a(t , pop , vlAdvancer , artDist , muHIV , ...
            kCD4 ,  maxRateM1 , maxRateM2 , maxRateF1 , maxRateF2 , disease , ...
            viral , gender , age , risk , k , hivInds , ...
            stepsPerYear , year) , tspan , pop(end , :));
        artTreatTracker(i , : , : , : , :  ,:) = artTreat;
        if any(pop(end , :) < 0)
            disp('After hiv')
            break
        end
    end
    
    [~ , pop , deaths(i , :) , vaxd(i , :)] = ode4xtra(@(t , pop) ...
        bornAgeDie(t , pop , ager , year , currStep , age , fertility , ...
        fertMat , fertMat2 , hivFertPosBirth ,hivFertNegBirth , hivFertPosBirth2 , ...
        hivFertNegBirth2 , deathMat , circMat , ...
        vaxerAger , vaxMat , MTCTRate , circStartYear , vaxStartYear ,...
        vaxRate , startYear , endYear, stepsPerYear) , tspan , pop(end , :));
    if any(pop(end , :) < 0)
        disp('After bornAgeDie')
        break
    end
    % add results to population vector
    popVec(i , :) = pop(end , :)';
end
popLast = popVec(end , :);

popVec = sparse(popVec); % compress population vectors

negSumLogL = likeFun(popVec , newCC , cinPos2014_obs , cinNeg2014_obs ,...
    hpv_hiv_2008_obs , hpv_hivNeg_2008_obs , hpv_hiv_obs , hpv_hivNeg_obs , ...
	hivPrevM_obs , hivPrevF_obs , disease , viral , gender , age , risk , ...
	hpvTypes , hpvStates , periods , startYear , stepsPerYear);


