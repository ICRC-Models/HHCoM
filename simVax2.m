% Accepts end year of analysis, number of tests, and population vector from
% calibrated model as inputs

% function simVax(lastYear , nTests, testParams , currPop)

close all; clear all; clc
lastYear = 2100;


disp('Start up')
% load population
popIn = load('H:\HHCoM_Results\kCC_Cin3_0.5Mult.mat');
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

w = ones(4 , 1) ./ 4;
kCC_Cin3_Orig = kCC_Cin3;
kCin2_Cin3_Orig = kCin2_Cin3;
kCin2_Cin1_Orig = kCin2_Cin1;
kCin1_Cin2_Orig = kCin1_Cin2;
kCin3_Cin2_Orig = kCin3_Cin2;
rNormal_Inf_Orig = rNormal_Inf;

for i = 1 : 3
    rNormal_Inf(: , i) = conv(rNormal_Inf_Orig(: , i) , w , 'same');
    rNormal_Inf(end - 1 : end , i) = conv(rNormal_Inf_Orig(end - 1 : end , i) , w , 'same');
    kCC_Cin3(: , i) = conv(kCC_Cin3_Orig(: , i) , w , 'same');
    kCC_Cin3(end - 1 : end , i) = kCC_Cin3_Orig(end - 1 : end , i);
%     kCin2_Cin3(: , i) = conv(kCin2_Cin3_Orig(: , i) , w , 'same');
%     kCin2_Cin3(end - 1 : end , i) = kCin3_Cin2_Orig(end - 1 : end , i);
    kCin3_Cin2(: , i) = conv(kCin3_Cin2_Orig(: , i) , w , 'same');
    kCin3_Cin2(end - 1 : end , i) = kCin3_Cin2_Orig(end - 1 : end , i);
    kCin1_Cin2(: , i) = conv(kCin1_Cin2_Orig(: , i) , w , 'same');
    kCin1_Cin2(end - 1 : end , i) = kCin1_Cin2_Orig(end - 1 : end , i);
    kCin2_Cin1(: , i) = conv(kCin2_Cin1_Orig(: , i) , w , 'same');
    kCin2_Cin1(end - 1 : end , i) = kCin2_Cin1_Orig(end - 1 : end , i);
end
% rNormal_Inf(: , 1) = rNormal_Inf(: , 1) .* 1.2;
% rNormal_Inf(: , 2) = rNormal_Inf(: , 2) .* 0.8;
% rNormal_Inf(: , 3) = rNormal_Inf(: , 3) .* 0.9;
rNormal_Inf = rNormal_Inf .* 0.85;
% c2c1Mults = 1.5 .* c2c1Mults;
% kCin1_Cin2(1 : end , :) = 1.2 * kCin1_Cin2(1 : end , :);
kCin2_Cin1(6 : end , :) = 1.25 * kCin2_Cin1(6 : end , :); 
kCin3_Cin2(10 : end , :) = 2.8 * kCin3_Cin2(10 : end , :);
kCin2_Cin3 = 0.5 .* kCin2_Cin3;
kCC_Cin3(7 : end , :) = 4 .* kCC_Cin3(7 : end , :);
muCC = min(muCC .* 12 , 0.99); % convert cervical cancer mortality rate from yearly to monthly
%     fImm(4 : age) = 1; % RR(0.75; 0.5 , 0.92) fraction fully protected by immunity based on RR of natural immunity (Beachler, 2017)
artHpvMult = hpv_hivMult(1 , :) * 0.25;
perPartnerHpv = 0.1; % high risk HPV transmission risk per month
rImmuneHiv = 3 ./ hpv_hivClear;
fImm(1 : age) = 1; % all infected individuals who clear HPV get natural immunity
lambdaMultImm(1 : 4) = 1 - 0.01;
lambdaMultImm(5 : 10) = 1 - logspace(log10(0.01) , log10(0.1) , 6);
lambdaMultImm(11 : age) = lambdaMultImm(10);
lambdaMultVax = ones(age , 2);
%%%%%

% load('HPV_calib.dat')
% kCin2_Cin3(: , 1) = HPV_calib(1 : age);
% kCin3_Cin2(: , 1) = HPV_calib(age + 1 : 2 * age);
% kCC_Cin3(: , 1) = HPV_calib(2 * age + 1 : 3 * age);
% kCin2_Cin3(: , 2) = HPV_calib(3 * age + 1 : 4 * age);
% kCin3_Cin2(: , 2) = HPV_calib(4 * age + 1 : 5 * age);
% kCC_Cin3(: , 2) = HPV_calib(5 * age + 1 : 6 * age);
% kCin2_Cin3(: , 3) = HPV_calib(6 * age + 1 : 7 * age);
% kCin3_Cin2(: , 3) = HPV_calib(7 * age + 1 : 8 * age);
% kCC_Cin3(: , 3) = HPV_calib(8 * age + 1 : 9 * age);
% rImmuneHiv = HPV_calib(9 * age + 1 : 9 * age + 1 + 3);
% c3c2Mults = HPV_calib(9 * age + 5 : 9 * age + 8);
% c2c1Mults = HPV_calib(9 * age + 9 : 9 * age + 12);
% artHpvMult = HPV_calib(9 * age + 13);
% perPartnerHpv = HPV_calib(9 * age + 14);

load('calibInitParams')
load('HPV_calib3.dat')
for i = 1 : 3
    kCin2_Cin3(: , i) = HPV_calib3(i) .* kCin2_Cin3(: , i);
    kCin3_Cin2(: , i) = HPV_calib3(3 + i) .* kCin3_Cin2(: , i);
    kCC_Cin3(: , i) = HPV_calib3(6 + i) .* kCC_Cin3(: , i);
end

rImmuneHiv = HPV_calib3(10 : 13);
c3c2Mults = HPV_calib3(14 : 17);
c2c1Mults = HPV_calib3(18 : 21);
artHpvMult = HPV_calib3(22);
perPartnerHpv= HPV_calib3(23);
lambdaMultImm = HPV_calib3(24 : 39);
%%%%%

c = fix(clock);
currYear = c(1); % get the current year
vaxEff = 0.8;
t_linearWane = 20; % pick a multiple of 5
k_wane = - vaxEff / t_linearWane;
lambdaMultVax_Arr = {zeros(age , 2) , zeros(age , 2) , zeros(age , 2) , ...
    zeros(age , 2)};
% 20 year waning
lambdaMultVax_Arr{1}(3 : 6 , 1) = vaxEff;
lambdaMultVax_Arr{1}(7 : 10 , 1) = vaxEff + k_wane * (5 : 5 : t_linearWane);
% 15 year waning
lambdaMultVax_Arr{2}(3 : 5 , 1) = vaxEff;
lambdaMultVax_Arr{2}(6 : 9 , 1) = vaxEff + k_wane * (5 : 5 : t_linearWane);
% 10 year waning
lambdaMultVax_Arr{3}(3 : 4 , 1) = vaxEff;
lambdaMultVax_Arr{3}(5 : 8 , 1) = vaxEff + k_wane * (5 : 5 : t_linearWane);
% 0 year waning
lambdaMultVax_Arr{4}(3 : age , 1) = vaxEff;

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
kCC_Cin3_Orig = kCC_Cin3;
kCC_Cin3 = kCC_Cin3_Orig .*0.5;
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
        currStep = round((s(i) + currYear) * stepsPerYear);
        tspan = [s(i) , s(i + 1)]; % evaluate diff eqs over one time interval
        popIn = popVec(i - 1 , :);
        if hpvOn
            hystOption = 'on';
            [~ , pop , newCC(i , : , : , : , :) , ccDeath(i , : , : , : , :)] ...
                = ode4xtra(@(t , pop) ...
                hpv(t , pop , immuneInds , infInds , cin1Inds , ...
                cin2Inds , cin3Inds , normalInds , ccInds , ccRegInds , ccDistInds , ...
                kInf_Cin1 , kCin1_Cin2 , kCin2_Cin3 , ...
                kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , kCin1_Inf  ,...
                rNormal_Inf , hpv_hivClear , c3c2Mults , ...
                c2c1Mults , fImm , kRL , kDR , muCC , disease , viral , age , hpvTypes , ...
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
            gar , perPartnerHpv , lambdaMultImm , lambdaMultVax , artHpvMult , epsA_vec , epsR_vec , yr , modelYr1 , ...
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
        
        
        [~ , pop , deaths(i , :)] = ode4xtra(@(t , pop) bornAgeDie(t , pop , ...
            ager , year , currStep , age , fertility , fertMat , hivFertPosBirth ,...
            hivFertNegBirth , deathMat , circMat , vaxerAger , MTCTRate , circStartYear , ...
            vaxStartYear , vaxRate , startYear , endYear , stepsPerYear) , tspan , pop(end , :));
        if any(pop(end , :) < 0)
            disp('After bornAgeDie')
            break
        end
        % add results to population vector
        popVec(i , :) = pop(end , :);
    end
    popLast = popVec(end , :);
    popVec = sparse(popVec); % compress population vectors
    filename = ['0.5 CIN3_CC' , 'Vax_' , num2str(vaxRate) , '_wane_' , ...
        num2str(k_wane) , '.mat']; %sprintf('test_output%d.mat' , n);
    parsave(filename , tVec ,  popVec , newHiv ,...
        newImmHpv , newVaxHpv , newHpv , deaths , hivDeaths , ccDeath , ...
        newCC , artTreatTracker , currYear , lastYear , popLast);
    
    %       savdir = 'H:\HHCoM Results';
    %     parsave(fullfile(savdir , 'resultsVax50') , 'tVec' ,  'popVec' , 'newHiv' ,...
    %         'newImmHpv' , 'newVaxHpv' , 'newHpv' , 'hivDeaths' , 'vaxRate' ,...
    %         'newCC' , 'artTreatTracker' , 'startYear' , 'endYear' , 'popLast');
end
disp('Done')
%%
simVaxResultOut()