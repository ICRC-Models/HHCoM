function [deaths, newCC, ccHealthState, symptomatic, treatment] = debugTreatment_sympCalibration_historicalSim(vaxResultInd , sceNum , fileNameNums, fileInds, deaths, ccHealthState, newCC, symptomatic, treatment , j)
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
    c3c2Mults , c2c1Mults , muCC , muCC_ud , muCC_d , kRL , kDR , artHpvMult , ...
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
    hystMult , ...
    hivNegNonVMMCinds , hivNegVMMCinds , vlAdvancer , ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , fertMat3 , hivFertPosBirth3 , hivFertNegBirth3 , ...
    fertMat4 , hivFertPosBirth4 , hivFertNegBirth4 , ...
    dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ...
    dFertPos3 , dFertNeg3  , dFertMat3, d_partnersMmult, riskAdj, d_riskAdj, ...
    deathMat , deathMat2 , deathMat3 , deathMat4 , deathMat5,...
    dDeathMat , dDeathMat2 , dDeathMat3 , dDeathMat4, dMue , ...
    ccLochpvVaxIndsFrom_treat , ...
    ccReghpvVaxInds_treat , ccDisthpvVaxInds_treat] = loadUp2(1 , 0 , [] , [] , []);

% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 1.5)

lastYear = 1950;  % ***SET ME***: last year of simulation (use 2122 for SA screening analysis)

%% Setting file names, initializing variables 
% nRuns = length(fileInds);

nRuns = 1; 

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
baseFileName = ['toNow_22Apr20Ph2V11_testingTreatment_changeSympProb_1925to1950_' sceNum '_1']; % ***SET ME***: name for simulation output file
% Looping length 
% loopSegments = {0 , round(nRuns/2) , nRuns};
% loopSegmentsLength = length(loopSegments);

%% Initialize screening inputs

% Treatment retention (proportion who return and comply with treatment)
% Note: ideally these parameters would be fed into the script via loadUp2.m
% These are not technically outputs of loadup2, but defined within loadup2,
% but not available outside of loadup2. 
% Retain is proportion not lost to follow up and return for the treatment. 
cryoRetain = 0.51; % with three-visit algorithm (cytology + colpo + cryotherapy treatment)
leepRetain = 0.80; % LEEP
thrmlRetain = 0.95; % thermal ablation
ccRetain = 0.40; % cancer treatment
eligLeep = [0.0 , 0.1 , 0.3]; % percent referred to/ eligible for LEEP (CIN1 , CIN2 , CIN3)

% if ((str2num(sceNum) == 0) || (str2num(sceNum) == 1)) % Scenarios 0 or 1
%     % Screening paper cytology algorithm
%     screenAlgs = spCyto;
%     screenAlgs.genTypBool = 0; % whether or not method uses HPV genotyping, only looks for vaccine types, only get treated if you have the high risk types
%     % proportion who return for treatment (susceptible/immune/infected/CIN1 , CIN2 , CIN3 , CC)
%     screenAlgs.leepRetain = zeros(1,4);
%     screenAlgs.cryoRetain = [0.0 , cryoRetain , cryoRetain , ccRetain];     
%     screenAlgs.thrmlRetain = zeros(1,4);  
% elseif ((str2num(sceNum) == 2) || (str2num(sceNum) == 3)) % Scenarios 2 or 3
%     % Screening paper HPV DNA -and-treat algorithm
%     screenAlgs = spHpvDna;
%     screenAlgs.genTypBool = 0;
%     % proportion who return for treatment (susceptible/immune/infected/CIN1 , CIN2 , CIN3 , CC)
%     screenAlgs.leepRetain = [[eligLeep.*leepRetain] , ccRetain]; 
%     screenAlgs.cryoRetain = zeros(1,4);
%     screenAlgs.thrmlRetain = [[(1-eligLeep).*thrmlRetain] , ccRetain]; 
% elseif ((str2num(sceNum) == 4) || (str2num(sceNum) == 5)) % Scenarios 4 or 5
%     % Screening paper HPV DNA+genotyping -and-treat algorithm
%     screenAlgs = spGentyp;
%     screenAlgs.genTypBool = 1;
%     % proportion who return for treatment (susceptible/immune/infected/CIN1 , CIN2 , CIN3 , CC)
%     screenAlgs.leepRetain = [[eligLeep.*leepRetain] , ccRetain]; 
%     screenAlgs.cryoRetain = zeros(1,4);
%     screenAlgs.thrmlRetain = [[(1-eligLeep).*thrmlRetain] , ccRetain]; 
% elseif ((str2num(sceNum) == 6) || (str2num(sceNum) == 7)) % Scenarios 6 or 7
%     % Screening paper AVE -and-treat algorithm
%     screenAlgs = spAve;
%     screenAlgs.genTypBool = 0;
%     % proportion who return for treatment (susceptible/immune/infected/CIN1 , CIN2 , CIN3 , CC)
%     screenAlgs.leepRetain = [[eligLeep.*leepRetain] , ccRetain]; 
%     screenAlgs.cryoRetain = zeros(1,4);
%     screenAlgs.thrmlRetain = [[(1-eligLeep).*thrmlRetain] , ccRetain]; 
% elseif ((str2num(sceNum) == 8) || (str2num(sceNum) == 9)) % Scenarios 8 or 9
%     % Screening paper HPV DNA + AVE triage -and-treat algorithm
%     screenAlgs = spHpvAve;
%     screenAlgs.genTypBool = 0;
%     % proportion who return for treatment (susceptible/immune/infected/CIN1 , CIN2 , CIN3 , CC)
%     screenAlgs.leepRetain = [[eligLeep.*leepRetain] , ccRetain]; 
%     screenAlgs.cryoRetain = zeros(1,4);
%     screenAlgs.thrmlRetain = [[(1-eligLeep).*thrmlRetain] , ccRetain]; 
% end

% Match the scenario to the ages that are screened for HIV pos and HIV neg 
if ismember(str2num(sceNum), [0 1 2 4 6 8])
    sceScreenAges = [8]; 
    sceScreenInd = {1}; 
    sceScreenIndAge = [8]; 
else
    sceScreenAges = [8 10 6 7 8 9 10]; % how screening ages are arranged in newScreen 
    sceScreenInd = {[1, 5] [2, 7] 3 4 6}; % these are the indices in newScreen. 2 and 7 should be summed. 1 and 5 should be summed. 
    sceScreenIndAge = [8 10 6 7 9]; 
end 

%% START FOR LOOP FOR RUNNING PARAMETERS 

% the for and parfor basically loops through the first half of the
% parameters, then the second half. but parfor allows them to run
% independently of each other. 
% for k = 1 : loopSegmentsLength-1
%     for j = loopSegments{k}+1 : loopSegments{k+1} % trying just a for loop for now

% k=1; % temporarily only running a single parameter to test 
j=1; 

% Load results
pathModifier = [baseFileName];
% nSims = size(dir([pwd , '/HHCoM_Results/' , pathModifier, '/' , '*.mat']) , 1);
curr = load([pwd , '/HHCoM_Results/' , pathModifier, '.mat']); % ***SET ME***: name for historical run output file
% vaxResult = cell(nSims , 1);
% resultFileName = [pwd , '/HHCoM_Results/' , pathModifier, '/' , 'vaxSimResult'];

% load results from vaccine run into cell array
% vaxResult{n} = load([resultFileName , num2str(n), '.mat']);

% concatenate vectors/matrices of population up to current year to population
% matrices for years past current year
% curr is historical model results 
% vaxResult is future model results
% this section of code combines the historical results with future
% results
% notice for vaxResult you start at row 2. likely because of
% 2021 being double counted in both. 
vaxResult{n}.popVec = [curr.popVec(1 : end  , :)]; % consolidating historical population numbers with future
vaxResult{n}.ccDeath = [curr.ccDeath(1 : end , : , : , :)]; % consolidating historical CC death #s with future... etc.
vaxResult{n}.ccDeath_treat = [curr.ccDeath_treat(1:end, :, :, :)];
vaxResult{n}.ccDeath_untreat = [curr.ccDeath_untreat(1:end, :, :, :)];
vaxResult{n}.ccDeath_treat_stage = [curr.ccDeath_treat_stage(1:end, :, :, :, :)];
vaxResult{n}.newCC = [curr.newCC(1 : end , : , : , :)];
vaxResult{n}.deaths = [curr.deaths(1 : end, 1)];
vaxResult{n}.ccSymp = [curr.ccSymp(1 : end, :, :, :, :, :, :, :)];
vaxResult{n}.ccTreat = [curr.ccTreat(1 : end, :, :, :, :, :, :, :)];
vaxResult{n}.newScreen = [curr.newScreen(1 : end, :, :, :, :, :, :, :)];

% DEATHS **************************************
    ccDeath = zeros(nTimepoints+1, age); 
    ccDeath_treat = zeros(nTimepoints+1, age); 
    ccDeath_untreat = zeros(nTimepoints+1, age); 
    ccDeath_treat_stage = zeros(nTimepoints+1, age, 6); 

    for a = 1 : age 
        for stage = 1 : 6
            ccDeath(:, a) = sum(sum(sum(vaxResult{n}.ccDeath(:, :, a, :),2),3),4); 
            ccDeath_treat(:, a) = sum(sum(sum(vaxResult{n}.ccDeath_treat(:, :, a, :),2),3),4); 
            ccDeath_untreat(:, a) = sum(sum(sum(vaxResult{n}.ccDeath_untreat(:, :, a, :),2),3),4); 
            ccDeath_treat_stage(:, a, stage) = sum(sum(sum(sum(vaxResult{n}.ccDeath_treat_stage(:, :, a, :, stage),2),3),4),5); 
        end 
    end 

    % combine all death data into 3D matrix
    deaths(:, 1:age, 1, 7, j) = ccDeath;
    deaths(:, 1:age, 2, 7, j) = ccDeath_treat; 
    deaths(:, 1:age, 3, 7, j) = ccDeath_untreat; 
    deaths(:, 1:age, 4, 1:6, j) = ccDeath_treat_stage; 

% CC HEALTH STATES ****************************

    for a = 1 : age 
	    for x = 1 : endpoints 
            cc = toInd(allcomb(1:disease, 1:viral, 6, 6, x, 1:intervens, 2, a, 1:risk)); 

		    ccHealthState(1:end, a, x, j) = sum(vaxResult{n}.popVec(:, cc), 2);
	    end 
    end 

% NEW CERVICAL CANCER CASES ******************************

    for a = 1 : age 
        for d = 1 : disease 
            for param = 1 : 2
                newCC(:, d, a, param, j) = sum(sum(sum(vaxResult{n}.newCC(:, d, a, param),2),3),4);
            end 
        end 
    end 

% SYMPTOMATIC CASES **************************************

    for x = 1 : 3
        for a = 1 : age
            for p = 1 : 3
                symptomatic(:, x, a, p, j) = sum(sum(sum(sum(sum(sum(sum(vaxResult{n}.ccSymp(:, :, :, :, x, :, a, p),2),3),4),5),6),7),8);
            end 
        end
    end 

% SCREENING ************************************************

%     for x = 1 : 3
%         for aInd = 1 : length(sceScreenInd)
%             a = sceScreenInd{aInd}; % the actual indices of screening age
%                 for p = 1 : 2
%                     screening(:, a, x, p, j) = sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(:, :, :, :, x, aInd, p),2),3),4),5),6),7);  
%                 end 
%         end 
%     end 

% TREATMENT ************************************************

    for x = 1 : 3
        for a = 1 : age
                for p = 1 : 3
                    treatment(:, x, a, p, j) = sum(sum(sum(sum(sum(sum(sum(vaxResult{n}.ccTreat(:, :, :, :, x, :, a, p),2),3),4),5),6),7),8); 
                end 
        end 
    end 

    

