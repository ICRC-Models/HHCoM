% Accepts end year of analysis, number of tests, and population vector from
% calibrated model as inputs

% function simVax(lastYear , nTests, testParams , currPop)

close all; clear all; clc
lastYear = 2080;


disp('Start up')
% load population
popIn = load('H:\HHCoM_Results\toNow');
currPop = popIn.popLast;

% load variables
paramDir = [pwd , '\Params\'];
load([paramDir, 'general'])
load([paramDir,'mixInfectIndices'])
load([paramDir,'vlAdvancer'])
load([paramDir,'fertMat'])
load([paramDir,'hivFertMats'])
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
load([paramDir ,'cost_weights'])
%%%%%%%
fImm(1 : age) = 1; % all infected individuals who clear HPV get natural immunity
%%%%%
%% use calibrated parameters
load([paramDir,'calibInitParams'])
load([paramDir,'HPV_calib12.dat'])
for i = 1 : 3
    kCin1_Inf(: , i) = HPV_calib12(i) .* kCin1_Inf(: , i);
    rNormal_Inf(: , i) = HPV_calib12(3 + i) .* rNormal_Inf(: , i);
    kCC_Cin3(: , i) = HPV_calib12(6 + i) .* kCC_Cin3(: , i);
    kCin3_Cin2(: , i) = HPV_calib12(49 + i) .* kCin3_Cin2(: , i); 
end

kCC_Cin3(: , 2) = kCC_Cin3(: , 3);
kCin3_Cin2(: , 3) = 1.5 .* kCin3_Cin2(: , 3);
kCC_Cin3(: , 2 : 3) = kCC_Cin3(: , 2 : 3) .* 1.25;

for i = 0 : 2
    maleActs(: , i + 1) = maleActs(: , i + 1) .* HPV_calib12(53 + i);
    femaleActs(: , i + 1) = femaleActs(: , i + 1) .* HPV_calib12(56 + i);
end

for a = 1 : age
    betaHIVF2M(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_F2M , maleActs(a , :)')); % HIV(-) males
    betaHIVM2F(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_M2F , femaleActs(a , :)')); % HIV(-) females
end
betaHIVM2F = permute(betaHIVM2F , [2 1 3]); % risk, age, vl
betaHIVF2M = permute(betaHIVF2M , [2 1 3]); % risk, age, vl
rImmuneHiv = HPV_calib12(10 : 13);
c3c2Mults = HPV_calib12(14 : 17);
c2c1Mults = HPV_calib12(18 : 21);
artHpvMult = HPV_calib12(22);
perPartnerHpv= HPV_calib12(23);
lambdaMultImm = HPV_calib12(24 : 39);
hpv_hivClear = HPV_calib12(40 : 43);
hpvClearMult = HPV_calib12(44 : 47);
perPartnerHpv_lr = HPV_calib12(48);%0.1;
perPartnerHpv_nonV = HPV_calib12(49); %0.1;

vaxMat = ager .* 0;
maxRateM_vec = [0.45 , 0.45];% maxRateM_arr{sim};
maxRateF_vec = [0.65 , 0.65];% maxRateF_arr{sim};

maxRateM1 = 1 - exp(-maxRateM_vec(1));
maxRateM2 = 1 - exp(-maxRateM_vec(2));
maxRateF1 = 1 - exp(-maxRateF_vec(1));
maxRateF2 = 1 - exp(-maxRateF_vec(2));
load([paramDir,'fertMat'])
load([paramDir,'hivFertMats'])
load([paramDir,'fertMat2'])
load([paramDir,'hivFertMats2'])
load([paramDir , 'settings']) 

artHpvMult = hpv_hivMult(1,1);
for h = 2 : 3
    hpv_hivMult(: , h) = hpv_hivMult(: , 1);
end

%%%%%

c = fix(clock);
currYear = c(1); % get the current year
vaxEff = 0.8;
t_linearWane = 20; % pick a multiple of 5
k_wane = - vaxEff / t_linearWane;
lambdaMultVax_Arr = {zeros(age , 2) , zeros(age , 2) , zeros(age , 2) , ...
    zeros(age , 2)};
for h = 1% : 2
    % 20 year waning
    lambdaMultVax_Arr{1}(3 : 6 , h) = vaxEff;
    lambdaMultVax_Arr{1}(7 : 10 , h) = (vaxEff + k_wane * (5 : 5 : t_linearWane))';
    % 15 year waning
    lambdaMultVax_Arr{2}(3 : 5 , h) = vaxEff;
    lambdaMultVax_Arr{2}(6 : 9 , h) = (vaxEff + k_wane * (5 : 5 : t_linearWane))';
    % 10 year waning
    lambdaMultVax_Arr{3}(3 : 4 , h) = vaxEff;
    lambdaMultVax_Arr{3}(5 : 8 , h) = (vaxEff + k_wane * (5 : 5 : t_linearWane))';
    % 0 year waning
    lambdaMultVax_Arr{4}(3 : age , h) = vaxEff;
end

k_wane_d = [20 , 15 , 10 , 0];
vaxCover = [0 , 0.5 , 0.7, 0.9];
testParams = allcomb(vaxCover , 1 : length(lambdaMultVax_Arr));
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
vaxStartYear = currYear;
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
    vaxerAgerArray{n} = vaxerAger;
end

% for males (future version?)
%     susMale = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 1 : periods , 1 , a , 1 : risk));
%     vaxdMale = toInd(allcomb(1 : disease , 1 : viral , 5 , 6 , 1 : periods , 1 , a , 1 : risk));
%     vaxer(vaxdMale , susMale) = V(2 , a);
%     vaxer(susMale , susMale) = -V(2 , a);
vaxMatArray = cell(3 , 1);
for n = 1 : size(testParams , 1)
    vaxMat = ager .* 0;
    vaxRate = testParams(n , 1);
    at = @(x , y) sort(prod(dim)*(y-1) + x);
    fromAge = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 1 , ...
        2 , 2 , 1 : risk));
    toAge = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 1 , ...
        2 , 3 , 1 : risk));
    toAgeVaxd = toInd(allcomb(1 : disease , 1 : viral , 1 , 9 , 1 , ...
        2 , 3 , 1 : risk));
    vaxMat(at(toAge , fromAge)) = (1 - vaxRate) * ager(at(toAge , fromAge));
    vaxMat(at(toAgeVaxd , fromAge)) = vaxRate * ager(at(toAge , fromAge)) ;
end
%%
parfor n = 1 : nTests
    vaxerAger = vaxerAgerArray{n};
    k_wane = k_wane_d(testParams(n , 2));
    lambdaMultVax = 1 - lambdaMultVax_Arr{testParams(n , 2)};
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
        currStep = round((s(i) + currYear - 1975) * stepsPerYear);
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
                rNormal_Inf , hpv_hivClear , c3c2Mults , hpvClearMult , ...
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
        popVec(i , :) = pop(end , :);
    end
    popLast = popVec(end , :);
    popVec = sparse(popVec); % compress population vectors
    filename = ['Vax_' , num2str(vaxRate) , '_wane_' , ...
        num2str(k_wane) , '.mat']; %sprintf('test_output%d.mat' , n);
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
%%
simVaxResultOut_AgeStand()
