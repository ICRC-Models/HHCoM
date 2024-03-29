function [hpvHealthState, hivHealthState, newCC, totalPerAge, totalPerAgeNoCC, onART] = vaxCEA_multSims_CISNET(vaxResultInd , sceNum , fileNameNums, fileInds, hpvHealthState, hivHealthState, newCC, totalPerAge, totalPerAgeNoCC, onART)
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
    intervens , gender , age , risk , hpvTypeGroups , dim , k , toInd , ...
    annlz , ...
    ageSexDebut , mInit , fInit , partnersM , partnersF , maleActs , ...
    femaleActs , riskDist , fertility , fertility2 , fertility3 , fertility4 , ...
    mue , mue2 , mue3 , mue4 , epsA_vec , epsR_vec , ...
    yr , ...
    hivOn , betaHIV_mod , muHIV , kCD4 , ...
    hpvOn , beta_hpvVax_mod , beta_hpvNonVax_mod , fImm , rImmune , ...
    kCin1_Inf , kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , rNormal_Inf , kInf_Cin1 , ...
    kCin1_Cin2 , kCin2_Cin3 , lambdaMultImm , hpv_hivClear , rImmuneHiv , ...
    c3c2Mults , c2c1Mults , c2c3Mults , c1c2Mults , muCC , muCC_ud , muCC_d , kRL , kDR , artHpvMult , ...
    hpv_hivMult , maleHpvClearMult , ...
    condUse , screenYrs , hpvScreenStartYear , ...
    artYr , maxRateM , maxRateF , ...
    artYr_vec , artM_vec , artF_vec , minLim , maxLim , ...
    circ_aVec , vmmcYr_vec , vmmc_vec , vmmcYr , vmmcRate , ...
    hivStartYear , circStartYear , circNatStartYear , vaxStartYear , ...
    baseline , who , spCyto , spHpvDna , spGentyp , spAve , spHpvAve , ...
    circProtect , condProtect , MTCTRate , hyst , ...
    OMEGA , ...
    ccInc2012_dObs , ccInc2018_dObs , cc_dist_dObs , cin3_dist_dObs , ...
    cin1_dist_dObs , hpv_dist_dObs , cinPos2002_dObs , cinNeg2002_dObs , ...
    cinPos2015_dObs , cinNeg2015_dObs , hpv_hiv_dObs , hpv_hivNeg_dObs , ...
    hpv_hivM2008_dObs , hpv_hivMNeg2008_dObs , hivPrevM_dObs , hivPrevF_dObs , ...
    popAgeDist_dObs , totPopSize_dObs , ...
    stageDist_1997_dObs , ...
    hivCurr , ...
    gar , hivSus , hpvVaxSus , hpvVaxImm , hpvNonVaxSus , hpvNonVaxImm , ...
    toHiv , vaxInds , nonVInds , hpvVaxInf , hpvNonVaxInf , ...
    hivInds , ...
    cin3hpvVaxIndsFrom , ccLochpvVaxIndsTo , ccLochpvVaxIndsFrom , ...
    ccReghpvVaxInds , ccDisthpvVaxInds , cin3hpvNonVaxIndsFrom , ...
    ccLochpvNonVaxIndsTo , ccLochpvNonVaxIndsFrom , ccReghpvNonVaxInds , ...
    ccDisthpvNonVaxInds , cin1hpvVaxInds , cin2hpvVaxInds , cin3hpvVaxInds , ...
    cin1hpvNonVaxInds , cin2hpvNonVaxInds , cin3hpvNonVaxInds , normalhpvVaxInds , ...
    immunehpvVaxInds , infhpvVaxInds , normalhpvNonVaxInds , immunehpvNonVaxInds , ...
    infhpvNonVaxInds , fromVaxNoScrnInds , fromVaxScrnInds , toNonVaxNoScrnInds , ...
    toNonVaxScrnInds , ageInd , riskInd , ...
    kSymp , hystMult , ...
    hivNegNonVMMCinds , hivNegVMMCinds , ...
    vlAdvancer , ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , fertMat3 , hivFertPosBirth3 , hivFertNegBirth3 , ...
    fertMat4 , hivFertPosBirth4 , hivFertNegBirth4 , ...
    dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ...
    dFertPos3 , dFertNeg3 , dFertMat3 , deathMat , deathMat2 , deathMat3 , deathMat4 , ...
    dDeathMat , dDeathMat2 , dDeathMat3 , dMue , ...
    ccLochpvVaxIndsFrom_treat , ...
    ccReghpvVaxInds_treat , ccDisthpvVaxInds_treat , vaxEff] = loadUp2_CISNET_S0(1 , 0 , [] , [] , [] , 1 , 1); % ***SET ME***

% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 1.5)

lastYear = 2124;  % ***SET ME***: last year of simulation (use 2122 for SA screening analysis)

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
diseaseVec_vax = {[1:2] [3:7] 8}; %indices for hiv neg, hiv pos no treatment, and hiv pos on art 
% Results directory
resultsDir = [pwd , '/HHCoM_Results/'];
fileKey = {'sim1' , 'sim0'};
fileKeyNums = fileNameNums;
n = vaxResultInd;
baseFileName = ['22Apr20Ph2V112v57BaseVax_spCytoScreen_shortName_noVMMC_noCond_noHiv_SA-S' , sceNum , '_']; % ***SET ME***: name for simulation output file
% Looping length 
loopSegments = {0 , round(nRuns/2) , nRuns};
loopSegmentsLength = length(loopSegments);

%% START FOR LOOP FOR RUNNING PARAMETERS 

% the for and parfor basically loops through the first half of the
% parameters, then the second half. but parfor allows them to run
% independently of each other. 
% for k = 1 : loopSegmentsLength-1
%     for j = loopSegments{k}+1 : loopSegments{k+1} % trying just a for loop for now

k=1; % temporarily only running a single parameter to test 
j=1; 

% Load results
pathModifier = [baseFileName , fileInds{j}];
nSims = size(dir([pwd , '/HHCoM_Results/' , pathModifier, '/' , '*.mat']) , 1);
curr = load([pwd , '/HHCoM_Results/toNow_22Apr20Ph2V11BaseVax_spCytoScreen_noVMMC_noCond_noHiv_' , fileInds{j}]); % ***SET ME***: name for historical run output file
vaxResult = cell(nSims , 1);
resultFileName = [pwd , '/HHCoM_Results/' , pathModifier, '/' , 'vaxSimResult'];

% load results from vaccine run into cell array
vaxResult{n} = load([resultFileName , num2str(n), '.mat']);

% concatenate vectors/matrices of population up to current year to population
% matrices for years past current year
% curr is historical model results 
% vaxResult is future model results
% this section of code combines the historical results with future
% results
% notice for vaxResult you start at row 2. likely because of
% 2021 being double counted in both. 
vaxResult{n}.popVec = [curr.popVec(1 : end  , :); vaxResult{n}.popVec(2 : end , :)]; % consolidating historical population numbers with future
vaxResult{n}.newCC = [curr.newCC(1 : end , : , : , :); vaxResult{n}.newCC(2 : end , : , : , :)];
vaxResult{n}.newHpvVax = [curr.newHpvVax(1 : end , : , : , : , : , :); vaxResult{n}.newHpvVax(2 : end , : , : , : , : , :)]; % infected with vaccine type HPV
vaxResult{n}.newImmHpvVax = [curr.newImmHpvVax(1 : end , : , : , : , : , :); vaxResult{n}.newImmHpvVax(2 : end , : , : , : , : , :)];
vaxResult{n}.newHpvNonVax = [curr.newHpvNonVax(1 : end , : , : , : , : , :); vaxResult{n}.newHpvNonVax(2 : end , : , : , : , : , :)];
vaxResult{n}.newImmHpvNonVax = [curr.newImmHpvNonVax(1 : end , : , : , : , : , :); vaxResult{n}.newImmHpvNonVax(2 : end , : , : , : , : , :)];
vaxResult{n}.artTreatTracker = [curr.artTreatTracker(1 : end , :  , : , : , : , :); vaxResult{n}.artTreatTracker(2 : end , : , : , : , : , :)];
vaxResult{n}.tVec = [curr.tVec(1 : end), vaxResult{n}.tVec(2 : end)];

% HPV HEALTH STATES ************************************
% To calculate HPV prevalence 
% CIN and CC also count towards HPV prevalence 

    for a = 1 : age
        for dInd = 1 : length(diseaseVec_vax)
            d = diseaseVec_vax{dInd}; 
            vaxInds1 = toInd(allcomb(d, 1:viral, 2:6, 2:6, 1:endpoints, 1:intervens, 2, a, 1:risk)); %women infected with both vax and nonvax types
            vaxInds2 = toInd(allcomb(d, 1:viral, [1 7], 2:6, 1:endpoints, 1:intervens, 2, a, 1:risk)); %women infected with only nonvax types
            vaxInds3 = toInd(allcomb(d, 1:viral, 2:6, [1 7], 1:endpoints, 1:intervens, 2, a, 1:risk)); %women infected with only vax types
            vaxInds = [vaxInds1; vaxInds2; vaxInds3]; %combine the three inds
            hpvHealthState(1:end, a, dInd, j) = sum(vaxResult{n}.popVec(:, vaxInds), 2);
	    end
    end


% HIV HEALTH STATES ************************************
% To calculate HIV prevalence 

    for a = 1 : age
        for dInd = 1 : length(diseaseVec_vax)
            d = diseaseVec_vax{dInd}; 
            vaxInds = toInd(allcomb(d, 1:viral, 1:hpvVaxStates, 1:hpvNonVaxStates, 1:endpoints, 1:intervens, 2, a, 1:risk)); 
            hivHealthState(1:end, a, dInd, j) = sum(vaxResult{n}.popVec(:, vaxInds), 2); 
        end
    end 

% NEW CERVICAL CANCER CASES ******************************
% To calculate CC incidence rate and CC case counts 

    for a = 1 : age 
        for dInd = 1 : length(diseaseVec_vax)
            d = diseaseVec_vax{dInd}; 
            newCC(1:end, a, dInd, j) = sum(sum(sum(vaxResult{n}.newCC(:, d, a, :),2),3),4);
        end 
    end 

% TOTAL NUMBER OF PEOPLE PER AGE GROUP AND HIV DISEASE STATE ********************
% To calculate crude prevalence 

    for a = 1 : age 
        for dInd = 1 : length(diseaseVec_vax)
            d = diseaseVec_vax{dInd}; 
            vaxInds = toInd(allcomb(d, 1:viral, 1:hpvVaxStates, 1:hpvNonVaxStates, 1:endpoints, 1:intervens, 2, a, 1:risk)); 
            totalPerAge(1:end, a, dInd, j) = sum(vaxResult{n}.popVec(:, vaxInds), 2); 
        end 
    end 

% TOTAL NUMBER OF PEOPLE WITHOUT CC ************************
% To calculate CC incidence (denominator is people without CC)

nonCC = [1 2 3 4 5 7]; % hpvvaxstates or nonvaxstates that are not CC

    for a = 1 : age
        for dInd = 1 : length(diseaseVec_vax)
            d = diseaseVec_vax{dInd}; 
            vaxInds = toInd(allcomb(d, 1:viral, nonCC, nonCC, 1:endpoints, 1:intervens, 2, a, 1:risk)); 
            totalPerAgeNoCC(1:end, a, dInd, j) = sum(vaxResult{n}.popVec(:, vaxInds), 2); 
        end 
    end 

% TOTAL NUMBER OF PEOPLE ON ART ****************************
% To calculate female ART coverage

    for a = 1 : age
        for dInd = 1 : length(diseaseVec_vax)
            d = diseaseVec_vax{dInd}; 
            onART(1:end, a, dInd, j) = sum(sum(sum(sum(sum(vaxResult{n}.artTreatTracker(:, d, 1:viral, 2, a, 1:risk), 2),3),4),5),6);
        end
    end

% end % for loop end 

% end % function end 

end