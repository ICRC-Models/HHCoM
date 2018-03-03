% Main
% Runs simulation over the time period and time step specified by the
% user.
close all; clear all; clc
profile clear;
% [~ , startYear , endYear , stepsPerYr , IntSteps] = Menu();
% disp('Done with input')
%%
% choose whether to model hysterectomy
hyst = 'off';
% choose whether to model HIV
hivOn = 1;
% Choose whether to model HPV
hpvOn = 0;
if hpvOn
    disp('HPV module activated.')
end

if hivOn
    disp('HIV module activated')
end
c = fix(clock);
currYear = c(1); % get the current year

% Get parameter values and load model

disp('Initializing. Standby...')
disp(' ');

startYear = 1975;
endYear = 2050; %currYear;
years = endYear - startYear;
save('settings' , 'years' , 'startYear' , 'endYear')
% Load parameters and constants for main
paramDir = [pwd , '\Params\'];
load([paramDir,'general'])
%% Initial population
load([paramDir,'popData'])
load([paramDir,'hpvData'])
% load('initPop')
% simulation
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
    initPop(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3) = 0.006 / 2 .* ...
        initPop_0(1 , 1 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3); % initial HIV infected male (4% prevalence)
    initPop(1 , 1 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3) = ...
        initPop_0(1 , 1 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3) .* (1 - 0.006 / 2); % moved to HIV infected
    initPop(3 , 2 , 1 , 1 , 1 , 2 , 4 : 6 , 2 : 3) = 0.006 / 2 .*...
        initPop_0(1 , 1 , 1 , 1 , 1 , 2 , 4 : 6 , 2 : 3); % initial HIV infected female (4% prevalence)
    initPop(1 , 1 , 1 , 1 , 1 , 2 , 4 : 6 , 2 : 3) = ...
        initPop_0(1 , 1 , 1 , 1 , 1 , 2 , 4 : 6 , 2 : 3) .* (1 - 0.006 / 2); % moved to HIV infected
    
    if hpvOn
        initPopHiv_0 = initPop;
        % HPV infected HIV+
        % females
        initPop(3 , 2 , 1 , 1 , 1 , 2 , 4 : 6 , 1 : 3) = 0.25 .* ...
            initPopHiv_0(3 , 2 , 1 , 1 , 1 , 2 , 4 : 6 , 1 : 3);
        
        % males
        initPop(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 1 : 3) = 0.25 .* ...
            initPopHiv_0(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 1 : 3);
        
        for h = 2 : 4
            % females
            initPop(3 , 2 , h , 1 , 1 , 2 , 4 : 6 , 1 : 3) = 0.75 / 3 .* ...
                initPopHiv_0(3 , 2 , 1 , 1 , 1 , 2 , 4 : 6 , 1 : 3);
            % males
            initPop(3 , 2 , h , 1 , 1 , 1 , 4 : 6 , 1 : 3) = 0.75 / 3 .* ...
                initPopHiv_0(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 1 : 3);
        end
    end
end
assert(~any(initPop(:) < 0) , 'Some compartments negative after seeding HIV infections.')

if hpvOn
    % initPop(1 , 1 , 2 : 4 , 1 , 1 , : , 4 : 9 , :) = 2; % initial HPV hr and lr infecteds (test)
    infected = initPop_0(1 , 1 , 1 , 1 , 1 , : , 6 : 9 , :) * 0.10; % try 10% intial HPV prevalence among age groups 6 - 9 (sexually active)
    initPop(1 , 1 , 1 , 1 , 1 , : , 6 : 9 , :) = ...
        initPop_0(1 , 1 , 1 , 1 , 1 , : , 6 : 9 , :) - infected;
    infected45 = initPop_0(1 , 1 , 1 , 1 , 1 , : , 4 : 5 , :) * 0.20; %try 20% initial HPV prevalence among age groups 4 - 5 (more sexually active)
    initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 5 , :) = ...
        initPop_0(1 , 1 , 1 , 1 , 1 , : , 4 : 5 , :) - infected45;
    % HPV 16/18
    initPop(1 , 1 , 2 , 1 , 1 , : , 6 : 9 , :) = 0.7 * infected;
    initPop(1 , 1 , 2 , 1 , 1 , : , 4 : 5 , :) = 0.7 * infected45;
    initPop(1 , 1 , 2 , 3 , 1 , : , 6 : 13 , :) = ...
        initPop_0(1 , 1 , 1 , 1 , 1 , : , 6 : 13 , :) .* 0.07 * 0.7;
    initPop(1 , 1 , 2 , 4 , 1 , : , 6 : 13 , :) = ...
        initPop_0(1 , 1 , 1 , 1 , 1 , : , 6 : 13 , :) .* 0.03 * 0.7;
    
    % 4v and oHR
    for h = 3 : 4
        initPop(1 , 1 , h , 1 , 1 , : , 6 : 9 , :) = infected ./ 3;
        initPop(1 , 1 , h , 1 , 1 , : , 4 : 5 , :) = infected45 ./ 3;
        initPop(1 , 1 , h , 3 , 1 , : , 6 : 13 , :) = ...
            initPop_0(1 , 1 , 1 , 1 , 1 , : , 6 : 13 , :) .* 0.07 * 0.3 / 2;
        initPop(1 , 1 , h , 4 , 1 , : , 6 : 13 , :) = ...
            initPop_0(1 , 1 , 1 , 1 , 1 , : , 6 : 13 , :) .* 0.03 * 0.3 / 2;
    end
    initPop(1 , 1 , 1 , 1 , 1 , : , 6 : 13 , :) = ...
        initPop_0(1 , 1 , 1 , 1 , 1 , : , 6 : 13 , :) * 0.9;
    initPop = max(initPop , 0);
end
assert(~any(initPop(:) < 0) , 'Some compartments negative after seeding HPV infections.')

% Intervention start years
circStartYear = 1990;
vaxStartYear = 2100;

%% Simulation
disp('Start up')
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
at = @(x , y) sort(prod(dim)*(y-1) + x);
k_wane = 0;
vaxRate = 0;
vaxerAger = ager;
fImm(1 : age) = 1; % all infected individuals who clear HPV get natural immunity
% profile on
disp(' ')
% Initialize vectors
timeStep = 1 / stepsPerYear;
years = endYear - startYear;
s = 1 : timeStep : years + 1; % stepSize and steps calculated in loadUp.m
artDistMat = zeros(size(prod(dim) , 20)); % initialize artDistMat to track artDist over past 20 time steps
%performance tracking
runtimes = zeros(size(s , 2) - 2 , 1);


%% use calibrated parameters
load([paramDir , 'calibInitParams'])
% load('HPV_calib3.dat')
% for i = 1 : 3
%     kCin2_Cin3(: , i) = HPV_calib3(i) .* kCin2_Cin3(: , i);
%     kCin3_Cin2(: , i) = HPV_calib3(3 + i) .* kCin3_Cin2(: , i);
%     kCC_Cin3(: , i) = HPV_calib3(6 + i) .* kCC_Cin3(: , i);
% end
% 
% rImmuneHiv = HPV_calib3(10 : 13);
% rImmuneHiv(1) = 1;
% c3c2Mults = max(1 , HPV_calib3(14 : 17));
% c2c1Mults = max(1 , HPV_calib3(18 : 21));
% artHpvMult = 1;%HPV_calib3(22);
% perPartnerHpv= HPV_calib3(23);
% lambdaMultImm = HPV_calib3(24 : 39);
% kCin1_Inf = 0.8 .* kCin1_Inf;
kCin2_Cin1 = 0.8 .* kCin2_Cin1;
rImmuneHiv = 3 ./ hpv_hivClear; 
perPartnerHpv = 0.1;%HPV_calib(9 * age + 14);
% maxRateM_vec = [0.25 , 0.35];
% maxRateF_vec = [0.35 , 0.5];

circProtect = 0.6;
dapProtect = zeros(age , 1);
dapProtect(4) = 0.27;
dapProtect(5 : age) = 0.5;
% maxRateM_arr = {[0.9 , 0.9] , [0.75 , 0.9] , [0.75 , 0.95] , [0.75 , 0.75]};
% maxRateF_arr = {[0.9 , 0.9] , [0.75, 0.9] , [0.75 , 0.95] , [0.75 , 0.75]};
testParams = [0 , 0.10 , 0.20];
dapAgerArray = cell(3 , 1);
for n = 1 : length(testParams)
    for a = 3 : 8
        dapAger = ager;
        dapCover = testParams(n);
        at = @(x , y) sort(prod(dim)*(y-1) + x);
        fromAge = toInd(allcomb(1 , 1 , 1 , 1 , 1 : periods , ...
            2 , a , 1 : risk));
        fromAgeDapd = toInd(allcomb(7 , 1 , 1 , 1 , 1 : periods , ...
            2 , a , 1 : risk));
        toAge = toInd(allcomb(1 , 1 , 1 , 1 , 1 : periods , ...
            2 , a + 1 , 1 : risk));
        toAgeDapd = toInd(allcomb(7 , 1 , 1 , 1 , 1 : periods , ...
            2 , a + 1 , 1 : risk));
        dapAger(at(toAge , fromAge)) = (1 - dapCover) * ager(at(toAge , fromAge));
        dapAger(at(toAge , fromAgeDapd)) = (1 - dapCover) * ager(at(toAge , fromAgeDapd));
        dapAger(at(toAgeDapd , fromAge)) = dapCover * ager(at(toAgeDapd , fromAge));
        dapAger(at(toAgeDapd , fromAgeDapd)) = dapCover * ager(at(toAgeDapd , fromAgeDapd));
        dapAgerArray{n} = dapAger;
    end
end
tits = {'noDap' ,'dap10' , 'dap20'};
% startYear = 2013;
endYear = 2050;
% dapAger = ager;
%%
disp(['Simulating period from ' num2str(startYear) ' to ' num2str(endYear) ...
    ' with ' num2str(stepsPerYear), ' steps per year.'])
disp(' ')
disp([num2str(length(testParams)) , ' simulation(s) running...'])
disp(' ')
% progressbar('Simulation Progress')
parfor sim = 1 : length(testParams)
    %%
% Gender and age specific max ART rates to test
    maxRateM_vec = [0.7 , 0.7];% maxRateM_arr{sim};
    maxRateF_vec = [0.7 , 0.7];% maxRateF_arr{sim};

    maxRateM1 = 1 - exp(-maxRateM_vec(1));
    maxRateM2 = 1 - exp(-maxRateM_vec(2));
    maxRateF1 = 1 - exp(-maxRateF_vec(1));
    maxRateF2 = 1 - exp(-maxRateF_vec(2));
    
    baseCirc = 1;
    dapAger = dapAgerArray{sim};
    % vectors to track specific pop changes
    %     artDistList = LinkedList();
    popVec = spalloc(years / timeStep , prod(dim) , 10 ^ 8);
    popIn = reshape(initPop , prod(dim) , 1); % initial population to "seed" model
    newHiv = zeros(length(s) - 1 , gender , age , risk);
    newHpv = zeros(length(s) - 1 , gender , disease , age , risk);
    newImmHpv = newHpv;
    newVaxHpv = newHpv;
    newCC = zeros(length(s) - 1 , disease , viral , hpvTypes , age);
    ccDeath = newCC;
    hivDeaths = zeros(length(s) - 1 , gender , age);
    deaths = popVec;
    artTreatTracker = zeros(length(s) - 1 , disease , viral , gender , age , risk);
    popVec(1 , :) = popIn;
    tVec = linspace(startYear , endYear , size(popVec , 1));
    k = cumprod([disease , viral , hpvTypes , hpvStates , periods , gender , age]);
    artDist = zeros(disease , viral , gender , age , risk); % initial distribution of inidividuals on ART = 0
    
    
    for i = 2 : length(s) - 1
        tic
        year = startYear + s(i) - 1;
%         disp(num2str(year))
        currStep = round(s(i) * stepsPerYear);
        %         disp(['current step = ' num2str(startYear + s(i) - 1) ' ('...
        %             num2str(length(s) - i) ' time steps remaining until year ' ...
        %             num2str(endYear) ')'])
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
            ode4xtra(@(t , pop) mixInfectDap(t , pop , currStep , ...
            gar , perPartnerHpv , lambdaMultImm , lambdaMultVax , artHpvMult , epsA_vec , epsR_vec , yr , modelYr1 , ...
            circProtect , condProtect , dapProtect , condUse , actsPer , partnersM , partnersF , ...
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
                kCD4 , maxRateM1 , maxRateM2 , maxRateF1 , maxRateF2 , ...
                disease , viral , gender , age , risk , k , hivInds , ...
                stepsPerYear , year) , tspan , pop(end , :));
            artTreatTracker(i , : , : , : , :  ,:) = artTreat;
            if any(pop(end , :) < 0)
                disp('After hiv')
                break
            end
            
%             [~ , pop , hivDeaths(i , : , :) , artTreat] =...
%                 ode4xtra(@(t , pop) hiv(t , pop , vlAdvancer , artDist , muHIV , ...
%                 kCD4 , disease , viral , gender , age , risk , k , hivInds , ...
%                 stepsPerYear , year) , tspan , pop(end , :));
%             artTreatTracker(i , : , : , : , :  ,:) = artTreat;
%             if any(pop(end , :) < 0)
%                 disp('After hiv')
%                 break
%             end
            %             [~ , artTreat] = ode4x(@(t , artDist) treatDist(t , popCopy(end , :) , year) , tspan , artDist);
            %             if size(artDistList) >= 20
            %                 artDistList.remove(); % remove earlier artDist matrix more than "20 time steps old"
            %             else
            %                 artDistList.add(artTreat);
            %             end
            %             artDist = calcDist(artDistList);
        end

%         [~ , pop , deaths(i , :)] = ode4xtra(@(t , pop) bornAgeDie(t , pop , ...
%             ager , year , currStep , age , fertility , fertMat , hivFertPosBirth ,...
%             hivFertNegBirth , deathMat , circMat , circAger , MTCTRate , circStartYear , ...
%             vaxStartYear , vaxRate , startYear , endYear , stepsPerYear) , tspan , pop(end , :));
        
        [~ , pop , deaths(i , :)] = ode4xtra(@(t , pop) bornAgeDieDap(t , pop ,...
            ager , year , currStep , age , fertility , fertMat , hivFertPosBirth ,...
            hivFertNegBirth , deathMat , baseCirc , circMat , ...
            dapAger , MTCTRate , circStartYear , vaxStartYear , ...
            vaxRate , startYear , endYear, currYear , stepsPerYear) , tspan , pop(end , :));
        
        if any(pop(end , :) < 0)
            disp('After bornAgeDie')
            break
        end
        % add results to population vector
        popVec(i , :) = pop(end , :)';
        %         runtimes(i) = toc;
        %         progressbar(i/(length(s) - 1))
    end
    popLast = popVec(end , :);
    %     disp(['Reached year ' num2str(endYear)])
    popVec = sparse(popVec); % compress population vectors
    %For local runs
    % savdir = 'C:\Users\nicktzr\Google Drive\ICRC\CISNET\Results';
    % save(fullfile(savdir , 'to2017') , 'tVec' ,  'popVec' , 'newHiv' ,...
    %     'newImmHpv' , 'newVaxHpv' , 'newHpv' , 'hivDeaths' , ...
    %     'deaths' , 'newCC' , 'artTreatTracker' , 'startYear' , 'endYear' , 'popLast');
    % For cluster runs
    filename = tits{sim};
    parsaveHIV(filename , tVec ,  popVec , newHiv ,...
        newImmHpv , newVaxHpv , newHpv , deaths , ...
        hivDeaths , ccDeath , newCC , artTreatTracker , startYear ,...
        endYear , popLast);
    %%
end
disp(' ')

disp('Simulation complete.')

% profile viewer
%%
dapResults()
