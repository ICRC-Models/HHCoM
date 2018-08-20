% Main
% Runs simulation over the time period and time step specified by the
% user.
close all; clear all; clc
profile clear;
% [~ , startYear , endYear , stepsPerYr , IntSteps] = Menu();
% disp('Done with input')
%%
% Get parameter values and load model
disp('Initializing. Standby...')
disp(' ');
paramDir = [pwd , '\Params\'];
%% Model Specs
% choose whether to model hysterectomy
hyst = 'off';
% choose whether to model HIV
hivOn = 1;
% Choose whether to model HPV
hpvOn = 1;
if hpvOn
    disp('HPV module activated.')
end

if hivOn
    disp('HIV module activated')
end
c = fix(clock);
currYear = c(1); % get the current year
startYear = 1980;
endYear = currYear;
years = endYear - startYear;

load([paramDir , 'popData'])
load([paramDir , 'General'])
%% use calibrated parameters
load([paramDir , 'calibratedParams'])
maxRateM_vec = [0.4 , 0.4];% maxRateM_arr{sim};
maxRateF_vec = [0.5 , 0.5];% maxRateF_arr{sim};

maxRateM1 = maxRateM_vec(1);
maxRateM2 = maxRateM_vec(2);
maxRateF1 = maxRateF_vec(1);
maxRateF2 = maxRateF_vec(2);

%% Initial Population 
mInit = popInit(: , 1);

fInit = popInit(: , 2);
% use male risk distribution for both genders
riskDistF = riskDistM;
% partnersF = partnersM;

% use male partners for both genders
MpopStruc = riskDistM;
FpopStruc = riskDistF;

mPop = zeros(age , risk);
fPop = mPop;

for i = 1 : age
    mPop(i , :) = MpopStruc(i, :).* mInit(i) ./ 1.25;
    fPop(i , :) = FpopStruc(i, :).* fInit(i) ./ 1.25;
end

dim = [disease , viral , hpvTypes , hpvStates , periods , gender , age ,risk];
initPop = zeros(dim);
initPop(1 , 1 , 1 , 1 , 1 , 1 , : , :) = mPop;
initPop(1 , 1 , 1 , 1 , 1 , 2 , : , :) = fPop;
initPop_0 = initPop;
if hivOn
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

            for h = 2
                % females
                initPop(3 , 2 , h , 1 , 1 , 2 , 4 : 6 , 1 : 3) = 0.7 .* ...
                    initPopHiv_0(3 , 2 , 1 , 1 , 1 , 2 , 4 : 6 , 1 : 3);
               % males
                initPop(3 , 2 , h , 1 , 1 , 1 , 4 : 6 , 1 : 3) = 0.7 .* ...
                    initPopHiv_0(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 1 : 3);
            end
        end
end
assert(~any(initPop(:) < 0) , 'Some compartments negative after seeding HIV infections.')

if hpvOn
    infected = initPop_0(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) * 0.20; % 20% intial HPV prevalence among age groups 4 - 9 (sexually active)
    initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) = ...
        initPop_0(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) - infected;

    % Omni-HPV type (transition rates weighted by estimated prevalence in population)
    initPop(1 , 1 , 2 , 1 , 1 , : , 4 : 9 , :) = infected;
end
assert(~any(initPop(:) < 0) , 'Some compartments negative after seeding HPV infections.')

% Intervention start years
circStartYear = 1990;
vaxStartYear = 2017;

%% Simulation
disp('Start up')
at = @(x , y) sort(prod(dim)*(y-1) + x);
k_wane = 0;
vaxRate = 0;
vaxerAger = ager;
fImm(1 : age) = 1; % all infected individuals who clear HPV get natural immunity

profile on
disp(' ')
% Initialize vectors
timeStep = 1 / stepsPerYear;
years = endYear - startYear;
s = 1 : timeStep : years + 1; % stepSize and steps calculated in loadUp.m
artDistMat = zeros(size(prod(dim) , 20)); % initialize artDistMat to track artDist over past 20 time steps
%performance tracking
runtimes = zeros(size(s , 2) - 2 , 1);
import java.util.LinkedList
artDistList = LinkedList();
popVec = spalloc(years / timeStep , prod(dim) , 10 ^ 8);
popIn = reshape(initPop , prod(dim) , 1); % initial population to "seed" model
newHiv = zeros(length(s) - 1 , gender , age , risk);
newHpv = zeros(length(s) - 1 , gender , disease , age , risk);
newImmHpv = newHpv;
newVaxHpv = newHpv;
newCC = zeros(length(s) - 1 , disease , hpvTypes , age);
ccDeath = newCC;
ccTreated = zeros(length(s) - 1 , disease , hpvTypes , age , 3); % 3 cancer stages: local, regional, distant
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
    currStep = round(s(i) * stepsPerYear);
    disp(['current step = ' num2str(startYear + s(i) - 1) ' ('...
        num2str(length(s) - i) ' time steps remaining until year ' ...
        num2str(endYear) ')'])
    tspan = [s(i) , s(i + 1)]; % evaluate diff eqs over one time interval
    popIn = popVec(i - 1 , :);

    if hpvOn
        hystOption = 'on';
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
    
    %         [~ , pop , deaths(i , :) , vaxd(i , :)] = ode4xtra(@(t , pop) ...
    %             bornAgeDie(t , pop , ager , year , currStep , age , fertility , ...
    %             fertMat , fertMat2 , hivFertPosBirth ,hivFertNegBirth , hivFertPosBirth2 , ...
    %             hivFertNegBirth2 , deathMat , circMat , ...
    %             vaxerAger , vaxMat , MTCTRate , circStartYear , vaxStartYear ,...
    %             vaxRate , startYear , endYear, stepsPerYear) , tspan , pop(end , :));
    %         if any(pop(end , :) < 0)
    %             disp('After bornAgeDie')
    %             break
    %         end
    [~ , pop , deaths(i , :)] = ode4xtra(@(t , pop) ...
        bornAgeDieRisk(t , pop , year , currStep ,...
        gender , age , risk , fertility , fertMat , fertMat2 , hivFertPosBirth ,...
        hivFertNegBirth , hivFertPosBirth2 , hivFertNegBirth2 , deathMat , circMat , ...
        MTCTRate , circStartYear , ageInd , riskInd , riskDist , startYear , ...
        endYear, stepsPerYear) , tspan , pop(end , :));
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

savdir = [pwd , '\HHCoM_Results\'];%'H:\HHCoM_Results';
save(fullfile(savdir , 'toNow') , 'tVec' ,  'popVec' , 'newHiv' ,...
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
showResults()
