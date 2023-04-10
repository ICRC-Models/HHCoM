function [vax, newCC] = processModelResults_IPVC2023(vaxResultInd , sceNum , fileNameNums, fileInds, vax, newCC)
% Description: This function links with the script
% loopingCeaOverScenarios.m. It takes in initialized result variables and
% places the results into 3D matrices. Looks at death counts,
% screening/treatment counts, HPV/CC health state counts, HIV health state
% counts, total counts per age at each time point, and vaccination counts. 
% Example: vaxCEA_multSims_SACEA_CH(1 , '0' , {'0'})

%% Load parameters and results
paramDir = [pwd , '\Params\'];

[stepsPerYear , timeStep , startYear , currYear , endYear , ...
    years , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , ...
    intervens , gender , age , risk , hpvTypeGroups , dim , k , toInd, annlz , ...
    ageSexDebut , mInit , fInit , partnersM , partnersF , partnersMmult, maleActs , ...
    femaleActs , riskDist , fertility , fertility2 , fertility3 , fertility4,...
    mue , mue2 , mue3 , mue4 , mue5, epsA_vec , epsR_vec , yr , ...
    hivOn , betaHIV_mod , hiv_hpvMult, muHIV , kCD4 , ...
    hpvOn , beta_hpvVax_mod , beta_hpvNonVax_mod , fImm , rImmune , ...
    kCin1_Inf , kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , rNormal_Inf , kInf_Cin1 , ...
    kCin1_Cin2 , kCin2_Cin3 , lambdaMultImm , hpv_hivClear , rImmuneHiv , ...
    c3c2Mults , c2c1Mults , muCC , kRL , kDR , artHpvMult , ...
    hpv_hivMult , maleHpvClearMult , ...
    condUse , screenYrs , hpvScreenStartYear , waning , ...
    artYr , maxRateM , maxRateF , ...
    artYr_vec , artM_vec , artF_vec , minLim , maxLim , ...
    circ_aVec , vmmcYr_vec , vmmc_vec , vmmcYr , vmmcRate , ...
    hivStartYear , circStartYear ,circNatStartYear , vaxStartYear , baseline , cisnet , who , whob , ...
    circProtect , condProtect , MTCTRate , hyst , OMEGA , ...
    ccInc2012_dObs , cc_dist_dObs , cin3_dist_dObs , cin1_dist_dObs , ...
    hpv_dist_dObs , cinPos2007_dObs , cin1_2010_dObs , cin2_2010_dObs, ...
    hpv_hiv_dObs , hpv_hivNeg_dObs , hpv_all_dObs , hpv_hiv2009_dObs , ...
    hivPrevF_dObs , hivPrevM_dObs , hivPrevAll_dObs, popAgeDist_dObs , totPopSize_dObs , hivCurr , ...
    gar , hivSus , hpvVaxSus , hpvVaxImm , hpvNonVaxSus , hpvNonVaxImm , ...
    toHiv , vaxInds , nonVInds , hpvVaxInf , hpvNonVaxInf , hivInds , ...
    cin3hpvVaxIndsFrom , ccLochpvVaxIndsTo , ccLochpvVaxIndsFrom , ...
    ccReghpvVaxInds , ccDisthpvVaxInds , cin3hpvNonVaxIndsFrom , ...
    ccLochpvNonVaxIndsTo , ccLochpvNonVaxIndsFrom , ccReghpvNonVaxInds , ...
    ccDisthpvNonVaxInds , cin1hpvVaxInds , cin2hpvVaxInds , cin3hpvVaxInds , ...
    cin1hpvNonVaxInds , cin2hpvNonVaxInds , cin3hpvNonVaxInds , normalhpvVaxInds , ...
    immunehpvVaxInds , infhpvVaxInds , normalhpvNonVaxInds , immunehpvNonVaxInds , ...
    infhpvNonVaxInds , ageInd , riskInd , ...
    hivNegNonVMMCinds , hivNegVMMCinds , vlAdvancer , ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , fertMat3 , hivFertPosBirth3 , hivFertNegBirth3 , ...
    fertMat4 , hivFertPosBirth4 , hivFertNegBirth4 , ...
    dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ...
    dFertPos3 , dFertNeg3  , dFertMat3, d_partnersMmult, riskAdj, d_riskAdj, ...
    deathMat , deathMat2 , deathMat3 , deathMat4 , deathMat5,...
    dDeathMat , dDeathMat2 , dDeathMat3 , dDeathMat4, dMue] = loadUp2(1 , 0 , [] , [] , []);

% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 1.5)

lastYear = 2071;  % ***SET ME***: last year of simulation (use 2122 for SA screening analysis)

%% Setting file names, initializing variables 
nRuns = length(fileInds);

% Initialize model output plots
% Timespans
monthlyTimespan = [startYear : timeStep : lastYear]; % list all the timespans in a vector
monthlyTimespan = monthlyTimespan(1 : end-1); % remove the very last date
monthlyTimespanScreen = [endYear : timeStep : lastYear]; % screening time span starts at 2021
monthlyTimespanScreen = monthlyTimespanScreen(1:end-1); 
nTimepoints = length(monthlyTimespan);
nTimepointsScreen = length(monthlyTimespanScreen); 
% Population outputs
diseaseVec_vax = {[1:2], 3, 4, 5, 6, 7, 8}; % HIV negative grouped together, and then all the HIV positive states 
% Results directory
resultsDir = [pwd , '/HHCoM_Results/'];
fileKey = {'sim1' , 'sim0'};
fileKeyNums = fileNameNums;
n = vaxResultInd;
% ***SET ME***
baseFileName = ['08Apr_stochMod_2019vaxIntro_nona_to2023_S' , sceNum ]; % ***SET ME***: name for simulation output file
% Looping length 
loopSegments = {0 , round(nRuns/2) , nRuns};
loopSegmentsLength = length(loopSegments);

%% START FOR LOOP FOR RUNNING PARAMETERS 

% the for and parfor basically loops through the first half of the
% parameters, then the second half. but parfor allows them to run
% independently of each other. 
for k = 1 : loopSegmentsLength-1
    for j = loopSegments{k}+1 : loopSegments{k+1} % trying just a for loop for now

% k=1; % temporarily only running a single parameter to test 
% j=1; 

% Load results
% pathModifier = [baseFileName , fileInds{j}];
pathModifier = [baseFileName]; 
nSims = size(dir([pwd , '/HHCoM_Results/' , pathModifier, '/' , '*.mat']) , 1);
% ***SET ME***
curr = load([pwd , '/HHCoM_Results/toNow_08Apr23_stochMod_2019vaxIntro_coverage16_nona_to2023_' , fileInds{j}]); % ***SET ME***: name for historical run output file
vaxResult = cell(nSims , 1);
resultFileName = [pwd , '/HHCoM_Results/' , pathModifier, '/' , 'vaxSimResult' , fileInds{j}];

% load results from vaccine run into cell array
vaxResult{n} = load([resultFileName , '.mat']);

% concatenate vectors/matrices of population up to current year to population
% matrices for years past current year
% curr is historical model results 
% vaxResult is future model results
% this section of code combines the historical results with future
% results
% notice for vaxResult you start at row 2. likely because of
% 2021 being double counted in both. 
vaxResult{n}.popVec = [curr.popVec(1 : end  , :); vaxResult{n}.popVec(2 : end , :)]; % consolidating historical population numbers with future
% vaxResult{n}.ccDeath = [curr.ccDeath(1 : end , : , : , :) ; vaxResult{n}.ccDeath(2 : end , : , : , :)]; % consolidating historical CC death #s with future... etc.
vaxResult{n}.newCC = [curr.newCC(1 : end , : , : , :); vaxResult{n}.newCC(2 : end , : , : , :)];
% vaxResult{n}.deaths = [curr.deaths(1 : end, 1); vaxResult{n}.deaths(2 : end, 1)];
% vaxResult{n}.newHpvVax = [curr.newHpvVax(1 : end , : , : , : , : , :); vaxResult{n}.newHpvVax(2 : end , : , : , : , : , :)]; % infected with vaccine type HPV
% vaxResult{n}.newImmHpvVax = [curr.newImmHpvVax(1 : end , : , : , : , : , :); vaxResult{n}.newImmHpvVax(2 : end , : , : , : , : , :)];
% vaxResult{n}.newHpvNonVax = [curr.newHpvNonVax(1 : end , : , : , : , : , :); vaxResult{n}.newHpvNonVax(2 : end , : , : , : , : , :)];
% vaxResult{n}.newImmHpvNonVax = [curr.newImmHpvNonVax(1 : end , : , : , : , : , :); vaxResult{n}.newImmHpvNonVax(2 : end , : , : , : , : , :)];
% vaxResult{n}.newScreen = [vaxResult{n}.newScreen(1 : end , : , : , : , : , : , :)]; %[curr.newScreen(1 : end , : , : , : , : , : , : ); vaxResult{n}.newScreen(2 : end , : , : , : , : , : , :)];
% vaxResult{n}.newHiv = [curr.newHiv(1 : end , : , : , : , : , : , :); vaxResult{n}.newHiv(2 : end , : , : , : , : , : , :)];
% vaxResult{n}.hivDeaths = [curr.hivDeaths(1 : end , : , : , :); vaxResult{n}.hivDeaths(2 : end , : , : , :)];
% vaxResult{n}.artTreatTracker = [curr.artTreatTracker(1 : end , :  , : , : , : , :); vaxResult{n}.artTreatTracker(2 : end , : , : , : , : , :)];
% vaxResult{n}.tVec = [curr.tVec(1 : end), vaxResult{n}.tVec(2 : end)];
vaxResult{n}.vaxdSchool = [curr.vaxdSchool(:, :); vaxResult{n}.vaxdSchool(2:end, :)]; 

% VACCINATIONS ********************************
    vaxTemplate = zeros(nTimepoints, 1); 

    vaxTemplate(:, 1) = vaxResult{n}.vaxdSchool; 
    
    vax(:, 17, j) = vaxTemplate; % input into the "all age" category 17

% NEW CERVICAL CANCER CASES ******************************

    % vaxResult{n}.newCC % (time, disease, age, hpvTypeGroups (2))

    for a = 1 : age 
        newCC(:, a, j) = sum(sum(sum(vaxResult{n}.newCC(:, :, a, :),2),3),4);
    end 

end % for loop end 

end % function end 