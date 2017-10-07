% Accepts end year of analysis, number of tests, and population vector from
% calibrated model as inputs

% function simVax(lastYear , nTests, testParams , currPop)

close all; clear all; clc
lastYear = 2100;


disp('Start up')
% load population
popIn = load('H:\HHCoM Results\results');
currPop = popIn.popLast;

% load variables
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
load('vaxInds')
%%%%%%%
load('calibParams')
fImm(1 : age) = 1; % all infected individuals who clear HPV get natural immunity
lambdaMultImm(1 : 4) = 1 - 0.01;
lambdaMultImm(5 : 10) = 1 - logspace(log10(0.01) , log10(0.1) , 6);
lambdaMultImm(11 : age) = lambdaMultImm(10);
lambdaMultVax = 1 - (0.9 * 0.8);

c = fix(clock);
currYear = c(1); % get the current year
vaxEff = 0.9;
k_wane = - vaxEff / 20;
lambdaMultImm_Arr = {zeros(age , 1) , zeros(age , 1) , zeros(age , 1)};
lambdaMultImm_Arr{1}(3 : 6) = vaxEff;
lambdaMultImm_Arr{1}(7 : 10) = vaxEff + k_wane * [5 : 5 : 20];
lambdaMultImm_Arr{2}(3 : 5) = vaxEff;
lambdaMultImm_Arr{2}(6 : 9) = vaxEff + k_wane * [5 : 5 : 20];
lambdaMultImm_Arr{3}(3 : 4) = vaxEff;
lambdaMultImm_Arr{3}(5 : 8) = vaxEff + k_wane * [5 : 5 : 20];



vaxCover = [0 , 0.5 , 0.7, 0.9];
testParams = allcomb(vaxCover , k_wane);
% titles = {'noVax' , 'vax50' , 'vax70' , 'vax90'};
nTests = size(testParams , 1);

%%%%%%%
dim = [disease , viral , hpvTypes , hpvStates , periods , gender , age ,risk];
% run analyses
stepsPerYear = 6;
c = fix(clock);
currYear = c(1); % get the current year
% Initialize vectors
timeStep = 1 / stepsPerYear;
years = lastYear - currYear;
s = 1 : timeStep : years + 1;
hivOn = 1;
hpvOn = 1;
% Intervention start years
circStartYear = 1990;
vaxStartYear = 2017;
disp(['Simulating period from ' num2str(currYear) ' to ' num2str(lastYear) ...
    ' with ' num2str(stepsPerYear), ' steps per year.'])
vaxerAgerArray = cell(3 , 1);
% for n = 1 : length(testParams)
%     vaxer = speye(length(currPop));
%     vaxRate = testParams(n);
%     at = @(x , y) sort(prod(dim)*(y-1) + x);
%     susFemale = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 1 : periods , 2 , 3 , 1 : risk));
%     vaxdFemale = toInd(allcomb(1 : disease , 1 : viral , 1 , 9 , 1 : periods , 2 , 3 , 1 : risk));
%     vaxer(at(vaxdFemale , susFemale)) = 1/5 .* vaxRate;
%     vaxer(at(susFemale , susFemale)) = 1 - vaxRate .* 1/5;
%     vaxerArray{n} = vaxer;
% end
for n = 1 : size(testParams , 1)
    vaxerAger = ager;
    vaxRate = testParams(n , 1);
    at = @(x , y) sort(prod(dim)*(y-1) + x);
    fromAge = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 1 : periods , ...
        2 , 2 , 1 : risk));
    toAge = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 1 : periods , ...
        2 , 3 , 1 : risk));
    toAgeVaxd = toInd(allcomb(1 : disease , 1 : viral , 1 , 9 , 1 : periods , ...
        2 , 3 , 1 : risk));
    vaxerAger(at(toAge , fromAge)) = (1 - vaxRate) * ager(at(toAge , fromAge));
    vaxerAger(at(toAgeVaxd , fromAge)) = vaxRate * ager(at(toAge , fromAge)) ;
%     vaxer(at(toAge , fromAge)) = 1 - vaxRate;
    vaxerAgerArray{n} = vaxerAger;
end

    % for males (future version?)
    %     susMale = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 1 : periods , 1 , a , 1 : risk));
    %     vaxdMale = toInd(allcomb(1 : disease , 1 : viral , 5 , 6 , 1 : periods , 1 , a , 1 : risk));
    %     vaxer(vaxdMale , susMale) = V(2 , a);
    %     vaxer(susMale , susMale) = -V(2 , a);

%%
parfor n = 1 : nTests
    vaxerAger = vaxerAgerArray{n};
    k_wane = testParams(n , 2);
    vaxRate = testParams(n , 1);
    popVec = spalloc(years / timeStep , prod(dim) , 10 ^ 8);
    popIn = currPop; % initial population to "seed" model
    newHiv = zeros(length(s) - 1 , gender , age , risk);
    newHpv = zeros(length(s) - 1 , gender , disease , age , risk);
    newImmHpv = newHpv;
    newVaxHpv = newHpv;
    newCC = zeros(length(s) - 1 , disease , viral , hpvTypes , age);
    ccDeath = newCC;
    hivDeaths = zeros(length(s) - 1 , age);
    deaths = zeros(size(popVec));
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
        if hivOn
            [~ , pop , newHiv(i , : , : , :)] = ...
                ode4xtra(@(t , pop) mixInfectHIV(t , pop , currStep , ...
                gar , hivSus , toHiv , mCurr , fCurr , ...
                mCurrArt , fCurrArt ,epsA_vec , epsR_vec , yr , modelYr1 , ...
                circProtect , condProtect , condUse , actsPer , partnersM , partnersF , ...
                betaHIVF2M , betaHIVM2F , disease , viral , gender , age , risk , ...
                hpvStates , hpvTypes , k , periods , stepsPerYear , year) , tspan , popIn);
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
            vaxStartYear , vaxRate,  startYear , lastYear , stepsPerYear) , tspan , pop(end , :));
%         if year >= vaxStartYear
%             pop = vaxer * pop(end , :)';
%             pop = pop';
%         end
        if any(pop(end , :) < 0)
            disp('After bornAgeDie')
            break
        end
        % add results to population vector
        popVec(i , :) = pop(end , :);

    end
    popLast = popVec(end , :);
    popVec = sparse(popVec); % compress population vectors
    filename = ['Vax_' , num2str(vaxRate) , '_wane_' , ...
        num2str(k_wane)]; %sprintf('test_output%d.mat' , n);
    parsave(filename , tVec ,  popVec , newHiv ,...
        newImmHpv , newVaxHpv , newHpv , hivDeaths , ccDeath , ...
        newCC , artTreatTracker , currYear , lastYear , popLast);

%       savdir = 'H:\HHCoM Results';
%     parsave(fullfile(savdir , 'resultsVax50') , 'tVec' ,  'popVec' , 'newHiv' ,...
%         'newImmHpv' , 'newVaxHpv' , 'newHpv' , 'hivDeaths' , 'vaxRate' ,...
%         'newCC' , 'artTreatTracker' , 'startYear' , 'endYear' , 'popLast');
end
disp('Done')
