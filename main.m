% Main
% Runs simulation over the time period and time step specified by the
% user.
close all; clear all; clc
profile clear;
% [~ , startYear , endYear , stepsPerYr , IntSteps] = Menu();
% disp('Done with input')
%%
% calcInds(stepsPerYear); % uncomment if recalculating indices
% Get parameter values and load model
prompt = {'End year' , ...
    'Model HIV? (1 = yes , 0 = no)' , ...
    'Model HPV? (1 = yes , 0 = no)' , 'Load new data? (1 = yes , 0 = no)' , ...
    'Calibrate? (1 = yes , 0 = no)'};
tit = 'Model Inputs';
inputs = inputdlg(prompt , tit);
calibrate = str2double(inputs{5});
loadNew = str2double(inputs{4});
c = fix(clock);
currYear = c(1); % get the current year
if isempty(inputs)
    return
end
disp('Initializing. Standby...')
disp(' ');
% Choose whether to model HPV
hpvOn = str2double(inputs{3});
if hpvOn
    disp('HPV module activated.')
end
% choose whether to model hysterectomy
hyst = 'off';
% choose whether to model HIV
hivOn = str2double(inputs{2});
if hivOn
    disp('HIV module activated')
end
disp(' ')
startYear = 1980;
endYear = str2double(inputs{1});
years = endYear - startYear;
% stepsPerYear = 4; % must reload data if changed
save('settings' , 'years' , 'startYear' , 'endYear')
if loadNew
    prompt = {'Enter steps per year'};
    input = inputdlg(prompt , tit);
    stepsPerYear = str2double(input{1});
    disp('Loading data. Standby...')
    loadUp(stepsPerYear);
else
    disp('Skipped parameter load up.')
    disp(' ')
end
% Load parameters and constants for main
load('general')
%% Initial population
load('popData')
if calibrate % calibrate model
%     disp('Pulling variables to be calibrated...')
%     file = 'PopData.xlsx';
%     maleActs = xlsread(file , 'Demographics' , 'AJ10 : AL25');
%     femaleActs = xlsread(file , 'Demographics' , 'AJ32 : AL47');
%     epsA = xlsread(file , 'Demographics' , 'AO2:AO4'); % [year] <1998 , 2003 , >2010 ; force of infection mixing
%     epsR = xlsread(file , 'Demographics' , 'AP2:AP4'); % [year] <1998 , 2003 , >2010
    disp('Running immunity calibration')
    fImm(1 : 3) = 1;
    fImm(4 : age) = 0.58; % (0.48; 0.27 , 0.69) fraction fully protected by immunity based on RR of natural immunity (Beachler, 2017)
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    nonlcon = [];
    x0 = [fImm(:)];
    lb = [0.27 * ones(length(fImm(:)) , 1)];
    ub = [0.69 * ones(length(fImm(:)) , 1)];
    options  = optimoptions('patternsearch' , 'UseParallel' , true , ...
        'UseCompletePoll' , true , 'PlotFcn' , ...
        {@psplotbestf , @psplotfuncount , @psplotmeshsize} , ...
        'AccelerateMesh' , true);
    disp('Calibration started.')
    calibParams = patternsearch(@calibrator , x0 , A , b , Aeq , beq , lb , ub , nonlcon , options);
    disp('Calibration complete.')
    save('calibParams' , 'calibParams')
else % simulation
    mInit = popInit(: , 1);
    MsumInit = sum(mInit);

    fInit = popInit(: , 2);
    FsumInit = sum(fInit);

    MpopStruc = riskDistM;
    FpopStruc = riskDistF;

    mPop = zeros(age , risk);
    fPop = mPop;

    for i = 1 : age
        mPop(i , :) = MpopStruc(i, :).* mInit(i) / 1.13;
        fPop(i , :) = FpopStruc(i, :).* fInit(i) / 1.13;
    end

    dim = [disease , viral , hpvTypes , hpvStates , periods , gender , age ,risk];
    initPop = zeros(dim);
    initPop(1 , 1 , 1 , 1 , 1 , 1 , : , :) = mPop;
    initPop(1 , 1 , 1 , 1 , 1 , 2 , : , :) = fPop;

    if hpvOn
        % initPop(1 , 1 , 2 : 4 , 1 , 1 , : , 4 : 9 , :) = 2; % initial HPV hr and lr infecteds (test)
        infected = initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) * 0.1; % try 10% intial HPV prevalence among age groups 4 - 9 (sexually active)
        initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) = ...
            initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) - infected;
        for h = 2
            initPop(1 , 1 , 2 , 1 , 1 , : , 4 : 9 , :) = infected;
        end
        initPop = max(initPop , 0);
    end
    assert(~any(initPop(:) < 0) , 'Some compartments negative after seeding HPV infections.')

    if hivOn
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
    end
    assert(~any(initPop(:) < 0) , 'Some compartments negative after seeding HIV infections.')

    % Intervention start years
    circStartYear = 1990;
    vaxStartYear = 2017;

    %% Simulation
    disp('Start up')
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
    load('calibParams') % calibrated fImm
    fImm = calibParams;
    OMEGA(1 : 3) = 0;
    OMEGA(4 : age) = logspace(-log(1 - 0.05) , - log(1 - 0.4) , age - 3);
    hivPositiveArtAll = sort(toInd(allcomb(10 , 6 , 1 : hpvStates , 1 : hpvTypes , ...
        1 : periods , 1 : gender  , 1 : age , 1 : risk)));
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
    newHpv = newHiv;
    newCC = zeros(length(s) - 1 , disease , viral , hpvTypes , age);
    hivDeaths = zeros(length(s) - 1 , age);
    deaths = popVec;
    artTreatTracker = zeros(length(s) - 1 , disease , viral , gender , age , risk);
    popVec(1 , :) = popIn;
    tVec = linspace(startYear , endYear , size(popVec , 1));
    k = cumprod([disease , viral , hpvTypes , hpvStates , periods , gender , age]);
    artDist = zeros(disease , viral , gender , age , risk); % initial distribution of inidividuals on ART = 0
    disp(['Simulating period from ' num2str(startYear) ' to ' num2str(endYear) ...
        ' with ' num2str(stepsPerYear), ' steps per year.'])
    disp(' ')
    disp('Simulation running...')
    disp(' ')
    try
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
                [~ , pop , newHpv(i , : , : , :)] = ...
                    ode4xtra(@(t , pop) mixInfectHPV(t , pop , currStep , ...
                    gar , epsA_vec , epsR_vec , yr , modelYr1 , ...
                    circProtect , condProtect , condUse , actsPer , partnersM , partnersF , ...
                    hpv_hivMult , hpvSus , toHpv , disease , viral , gender , age , risk , hpvStates , hpvTypes , ...
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

%                 [~ , pop] = ode4x(@(t , pop) hpvTreat(t , pop , disease , viral , hpvTypes , age , ...
%                     periods , detCC , hivCC , muCC , ccRInds , ccSusInds , ...
%                     hystPopInds , screenFreq , screenCover , hpvSens , ccTreat , ...
%                     cytoSens , cin1Inds , cin2Inds , cin3Inds ,  normalInds , getHystPopInds ,...
%                     OMEGA , leep , hystOption , year) , tspan , pop(end , :));
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
            runtimes(i) = toc;
            progressbar(i/(length(s) - 1))
        end

        disp(['Reached year ' num2str(endYear)])
        popVec = sparse(popVec); % compress population vectors
        savdir = '/c/Users/nicktzr/Google Drive/ICRC/CISNET/Results';
	      save(fullfile(savdir , 'results') , 'results' , 'tVec' ,  'popVec' , 'newHiv' , 'newHpv' , 'hivDeaths' , ...
            'deaths' , 'newCC' , 'artTreatTracker' , 'startYear' , 'endYear');
        disp(' ')
        disp('Simulation complete.')
    catch
        disp('Simulation aborted.')
        error(['Error after ' num2str(startYear + s(i) - 1) ' ('...
            num2str(length(s) - i) ' time steps remaining until year ' ...
            num2str(endYear) ')'])
        save(fullfile(savdir , 'errorWorkspace'), 'errorWorkspace')
        disp('All variables saved to errorWorkspace for debugging.')
        return
    end
    profile viewer
    %%
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
    % %% load target values
    % load('targetVals')
    %% Show results
    showResults()
    %% Convert results to CSV
    % resultOut()
end
