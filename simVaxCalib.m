% Accepts end year of analysis, number of tests, and population vector from
% calibrated model as inputs

% function simVax(lastYear , nTests, testParams , currPop)

close all; clear all; clc
lastYear = 2081;

disp('Start up')

% load variables
paramDir = [pwd , '\Params\'];
load([paramDir, 'general'])

% load population
popIn = load('H:\HHCoM_Results\toNow');
currPop = popIn.popLast;
%%%%%%%
fImm(1 : age) = 1; % all infected individuals who clear HPV get natural immunity
%%%%%
%% use calibrated parameters
load([paramDir , 'General'])
load([paramDir , 'calibratedParams'])
maxRateM_vec = [0.45 , 0.45];% maxRateM_arr{sim};
maxRateF_vec = [0.6 , 0.6];% maxRateF_arr{sim};

maxRateM1 = 1 - exp(-maxRateM_vec(1));
maxRateM2 = 1 - exp(-maxRateM_vec(2));
maxRateF1 = 1 - exp(-maxRateF_vec(1));
maxRateF2 = 1 - exp(-maxRateF_vec(2));
%%%%%

% 90% efficacy against 70% of CC types, 100% efficacy against 70% of types, ...
% 100% efficacy against 90% of types
vaxEff = [0.9 * 0.7 , 0.7 , 0.9]; 
t_linearWane = 20; % pick a multiple of 5
k_wane = - vaxEff / t_linearWane;
vaxCover = [0 , 0.7 , 0.9];
testParams = allcomb(vaxCover , vaxEff);
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
vaxStartYear = currYear;

disp(['Simulating period from ' num2str(currYear) ' to ' num2str(lastYear) ...
    ' with ' num2str(stepsPerYear), ' steps per year.'])
fromNonV = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 1 , ...
    2 , 3 , 1 : risk));
toV = toInd(allcomb(1 : disease , 1 : viral , 1 , 9 , 1 , ...
    2 , 3 , 1 : risk));
import java.util.LinkedList
artDistList = LinkedList();
%%
parfor n = 1 : nTests
    vaxerAger = ager;% vaxerAgerArray{n};
    vaxMat = vaxMatArray{n};
    vaxEff = testParams(n , 2);
    lambdaMultVax = 1 - lambdaMultVaxMat(: , n);
    vaxRate = testParams(n , 1);
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
            [~ , pop , newCC(i , : , : , :) , ccDeath(i , : , : , :) , ...
                ccTreated(i , : , : , : , :)] ...
                = ode4xtra(@(t , pop) ...
                hpv(t , pop , immuneInds , infInds , cin1Inds , ...
                cin2Inds , cin3Inds , normalInds , ccInds , ccRegInds , ccDistInds , ...
                ccTreatedInds , kInf_Cin1 , kCin1_Cin2 , kCin2_Cin3 , ...
                kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , kCin1_Inf  ,...
                rNormal_Inf , hpv_hivClear , c3c2Mults , ...
                c2c1Mults , fImm , kRL , kDR , muCC , kCCDet , ...
                disease , viral , age , hpvTypes , ...
                rImmuneHiv , vaccinated , hystOption) , tspan , popIn);
            popIn = pop(end , :);
            if any(pop(end , :) <  0)
                disp('After hpv')
                break
            end
            
            %                 [~ , pop] = ode4x(@(t , pop) hpvTreat(t , pop , disease , viral , hpvTypes , age , ...
            %                     periods , detCC , hivCC , muCC , ccRInds , ccSusInds , ...
            %                     hystPopInds , screenFreq , screenCover , hpvSens , ccTreat , ...
            %                     cytoSens , cin1Inds , cin2Inds , cin3Inds ,  normalInds , getHystPopInds ,...
            %                     OMEGA , leep , hystOption , year) , tspan , pop(end , :));
        end
        
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
        
        if hivOn
            [~ , pop , hivDeaths(i , : , :) , artTreat] =...
                ode4xtra(@(t , pop) hiv2a(t , pop , vlAdvancer , artDist , muHIV , ...
                kCD4 ,  maxRateM1 , maxRateM2 , maxRateF1 , maxRateF2 , disease , ...
                viral , gender , age , risk , k , hivInds , ...
                stepsPerYear , year) , tspan , pop(end , :));
            artTreatTracker(i , : , : , : , :  ,:) = artTreat;
            artDistList.add(artTreat);
            if size(artDist) >= stepsPerYear * 5
                artDistList.remove(); % remove CD4 and VL distribution info for people initiating ART more than 5 years ago
            end
            artDist = calcDist(artDistList , disease , viral , gender , age , ...
                risk);
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
        
        
        [~ , pop , deaths(i , :)] = ode4xtra(@(t , pop) ...
            bornAgeDieRisk(t , pop , ager , year , currStep , gender , age , risk , fertility , ...
            fertMat , fertMat2 , hivFertPosBirth ,hivFertNegBirth , hivFertPosBirth2 , ...
            hivFertNegBirth2 , deathMat , circMat , ...
            vaxerAger , vaxMat , MTCTRate , circStartYear , ageInd ,...
            riskInd , riskDist , startYear , endYear, stepsPerYear) , tspan , pop(end , :));
        if any(pop(end , :) < 0)
            disp('After bornAgeDieRisk')
            break
        end
        
        fracVaxd = sum(pop(end , toV) , 2) / (sum(pop(end , fromNonV) , 2) + sum(pop(end , toV) , 2));
        if fracVaxd < vaxRate
            vaxCover = max(0 , (vaxRate - fracVaxd) ./ (1 - fracVaxd));
            pop(end , toV) = pop(end , toV) + pop(end , fromNonV) .* vaxCover;
            pop(end , fromNonV) = pop(end , fromNonV) .* (1 - vaxCover);
            vaxd(i , :) = sumall(pop(end , fromNonV).* vaxCover);
        end

        % add results to population vector
        popVec(i , :) = pop(end , :);
    end
    popLast = popVec(end , :);
    popVec = sparse(popVec); % compress population vectors
    filename = ['CEA_VaxCover_' , num2str(vaxRate) , '_Eff_' , ...
        num2str(vaxEff) , '.mat']; %sprintf('test_output%d.mat' , n);
    parsave(filename , tVec ,  popVec , newHiv ,...
        newImmHpv , newVaxHpv , newHpv , deaths , hivDeaths , ccDeath , ...
        newCC , artTreatTracker , vaxd , ccTreated , ...
        currYear , lastYear , popLast);
    
    %       savdir = 'H:\HHCoM Results';
    %     parsave(fullfile(savdir , 'resultsVax50') , 'tVec' ,  'popVec' , 'newHiv' ,...
    %         'newImmHpv' , 'newVaxHpv' , 'newHpv' , 'hivDeaths' , 'vaxRate' ,...
    %         'newCC' , 'artTreatTracker' , 'startYear' , 'endYear' , 'popLast');
end
disp('Done')
