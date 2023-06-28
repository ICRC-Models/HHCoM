function [vax, deaths, ccHealthState, hpvHealthState, newCC, totalPerAge, screenTreat, screenSympCCTreat] = vaxCEA_multSims_processResults_sce10(vaxResultInd , sceNum , fileNameNums, fileInds, vax, deaths, ccHealthState, hpvHealthState, newCC, totalPerAge, screenTreat, screenSympCCTreat)
% Description: This function links with the script
% loopingCeaOverScenarios.m. It takes in initialized result variables and
% places the results into 3D matrices. Looks at death counts,
% screening/treatment counts, HPV/CC health state counts, HIV health state
% counts, total counts per age at each time point, and vaccination counts. 

% This is for the Kenya 1-dose modeling study. 

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
    kSymp , hystMult , ...
    hivNegNonVMMCinds , hivNegVMMCinds , vlAdvancer , ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , fertMat3 , hivFertPosBirth3 , hivFertNegBirth3 , ...
    fertMat4 , hivFertPosBirth4 , hivFertNegBirth4 , ...
    dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ...
    dFertPos3 , dFertNeg3  , dFertMat3, d_partnersMmult, riskAdj, d_riskAdj, ...
    deathMat , deathMat2 , deathMat3 , deathMat4 , deathMat5,...
    dDeathMat , dDeathMat2 , dDeathMat3 , dDeathMat4, dMue , ...
    ccLochpvVaxIndsFrom_treat , ...
    ccReghpvVaxInds_treat , ccDisthpvVaxInds_treat , vaxEff] = loadUp2(1 , 0 , [] , [] , []);

lastYear = 2123; % manually set in futureSim

%% Setting file names, initializing variables 
nRuns = length(fileInds);

% Initialize model output plots
% Timespans
monthlyTimespan = [startYear : timeStep : lastYear]; % list all the timespans in a vector
monthlyTimespan = monthlyTimespan(1 : end-1); % remove the very last date
monthlyTimespanScreen = [2000 : timeStep : lastYear]; % screening time span starts at 2021
monthlyTimespanScreen = monthlyTimespanScreen(1:end-1); 
nTimepoints = length(monthlyTimespan);
nTimepointsScreen = length(monthlyTimespanScreen); 
% Population outputs
% diseaseVec_vax = {[1:2], 3, 4, 5, 6, 7, 8}; % HIV negative grouped together, and then all the HIV positive states 
% Results directory
resultsDir = [pwd , '/HHCoM_Results/'];
fileKey = {'sim1' , 'sim0'};
fileKeyNums = fileNameNums;
n = vaxResultInd;
baseFileName = ['VaccineKenya1DoseCea_S' , sceNum]; % ***SET ME***: name for simulation output file
% Looping length 
loopSegments = {0 , round(nRuns/2) , nRuns};
loopSegmentsLength = length(loopSegments);
fivYrAgeGrpsOn = 1;

% SCREENING
screenAlgorithm = 1; % ***SET ME***: screening algorithm to use (1 for baseline, 2 for CISNET, 3 for WHOa, 4 for WHOb)
hivPosScreen = 1; % ***SET ME***: 0 applies same screening algorithm (screenAlgorithm) for all HIV states; 1 applies screenAlgorithm to HIV+ and screenAlgorithmNeg to HIV-
screenAlgorithmNeg = 1; % ***SET ME***: If hivPosScreen=1, screening algorithm to use for HIV- persons (1 for baseline, 2 for CISNET, 3 for WHOa, 4 for WHOb) 

%% Screening
if (screenAlgorithm == 1)
    % Baseline screening algorithm
    screenAlgs{1} = baseline;
elseif (screenAlgorithm == 2)
    % CISNET screening algorithm
    screenAlgs{1} = cisnet;
elseif (screenAlgorithm == 3)
    % WHO screening algorithm - version a
    screenAlgs{1} = who;
elseif (screenAlgorithm == 4)
    % WHO screening algorithm - version b
    screenAlgs{1} = whob;
end

if hivPosScreen
    if (screenAlgorithmNeg == 1)
        % Baseline screening algorithm
        screenAlgs{2} = baseline;
    elseif (screenAlgorithmNeg == 2)
        % CISNET screening algorithm
        screenAlgs{2} = cisnet;
    elseif (screenAlgorithmNeg == 3)
        % WHO screening algorithm - version a
        screenAlgs{2} = who;
    elseif (screenAlgorithmNeg == 4)
        % WHO screening algorithm - version b
        screenAlgs{2} = whob;
    end
    screenAlgs{2}.screenCover_vec = cell(size(screenYrs , 1) - 1, 1); % save data over time interval in a cell array
    for i = 1 : size(screenYrs , 1) - 1          % interpolate dnaTestCover values at steps within period
        period = [screenYrs(i) , screenYrs(i + 1)];
        screenAlgs{2}.screenCover_vec{i} = interp1(period , screenAlgs{2}.screenCover(i : i + 1 , 1) , ...
            screenYrs(i) : timeStep : screenYrs(i + 1));
    end
    screenAlgs{1}.diseaseInds = [3 : 8];
    screenAlgs{2}.diseaseInds = [1 : 2];
else
    screenAlgs{1}.diseaseInds = [1 : disease];
end

screenAlgs{1}.screenCover_vec = cell(size(screenYrs , 1) - 1, 1); % save data over time interval in a cell array
for i = 1 : size(screenYrs , 1) - 1          % interpolate dnaTestCover values at steps within period
    period = [screenYrs(i) , screenYrs(i + 1)];
    screenAlgs{1}.screenCover_vec{i} = interp1(period , screenAlgs{1}.screenCover(i : i + 1 , 1) , ...
        screenYrs(i) : timeStep : screenYrs(i + 1));
end

%% Initialize screening inputs

% Treatment retention (proportion who return and comply with treatment)
% Note: ideally these parameters would be fed into the script via loadUp2.m
% These are not technically outputs of loadup2, but defined within loadup2,
% but not available outside of loadup2. 
% Retain is proportion not lost to follow up and return for the treatment. 
colpoRetain = 0.72;
cinTreatRetain = 0.5;
ccTreatRetain = 0.4;

screenAge = [35/max(1 , fivYrAgeGrpsOn*5)+1];

% Match the scenario to the ages that are screened for HIV pos and HIV neg 
% if ismember(str2num(sceNum), [0 1 2 4 6 8])
%     sceScreenAges = [8]; 
%     sceScreenInd = {1}; 
%     sceScreenIndAge = [8]; 
% else
%     sceScreenAges = [8 10 6 7 8 9 10]; % how screening ages are arranged in newScreen 
%     sceScreenInd = {[1, 5] [2, 7] 3 4 6}; % these are the indices in newScreen. 2 and 7 should be summed. 1 and 5 should be summed. 
%     sceScreenIndAge = [8 10 6 7 9]; 
% end 

%% START FOR LOOP FOR RUNNING PARAMETERS 

% the for and parfor basically loops through the first half of the
% parameters, then the second half. but parfor allows them to run
% independently of each other. 
for k = 1 : loopSegmentsLength-1
    for j = loopSegments{k}+1 : loopSegments{k+1} % trying just a for loop for now

% k=1; % temporarily only running a single parameter to test 
% j=1; 

% Load results
pathModifier = [baseFileName , fileInds{j}];
nSims = size(dir([pwd , '/HHCoM_Results/' , pathModifier, '/' , '*.mat']) , 1);
curr = load([pwd , '/HHCoM_Results/toNow_07Jun23_stochMod_baseline_2dose_nowane' , fileInds{j}]); % ***SET ME***: name for historical run output file
vaxResult = cell(nSims , 1);
resultFileName = [pwd , '/HHCoM_Results/' , baseFileName, '/' , 'vaxWaneSimResult']; % ***SET ME***: change add "Wane" if you are processing waning results

% load results from vaccine run into cell array
vaxResult{n} = load([resultFileName , fileInds{j}, '.mat']);

% concatenate vectors/matrices of population up to current year to population
% matrices for years past current year
% curr is historical model results 
% vaxResult is future model results
% this section of code combines the historical results with future
% results
% notice for vaxResult you start at row 2. likely because of
% 2023 being double counted in both. 
vaxResult{n}.popVec = [curr.popVec(1 : end  , :); vaxResult{n}.popVec(2 : end , :)]; % consolidating historical population numbers with future
vaxResult{n}.ccDeath_treat = [curr.ccDeath_treat(1 : end , : , : , :) ; vaxResult{n}.ccDeath_treat(2 : end , : , : , :)]; % consolidating historical CC death #s with future... etc.
vaxResult{n}.ccDeath_untreat = [curr.ccDeath_untreat(1 : end , : , : , :) ; vaxResult{n}.ccDeath_untreat(2 : end , : , : , :)]; 
vaxResult{n}.newCC = [curr.newCC(1 : end , : , : , :); vaxResult{n}.newCC(2 : end , : , : , :)];
vaxResult{n}.deaths = [curr.deaths(1 : end, 1); vaxResult{n}.deaths(2 : end, 1)];
vaxResult{n}.newHpvVax = [curr.newHpvVax(1 : end , : , : , : , : , :); vaxResult{n}.newHpvVax(2 : end , : , : , : , : , :)]; % infected with vaccine type HPV
vaxResult{n}.newImmHpvVax = [curr.newImmHpvVax(1 : end , : , : , : , : , :); vaxResult{n}.newImmHpvVax(2 : end , : , : , : , : , :)];
vaxResult{n}.newHpvNonVax = [curr.newHpvNonVax(1 : end , : , : , : , : , :); vaxResult{n}.newHpvNonVax(2 : end , : , : , : , : , :)];
vaxResult{n}.newImmHpvNonVax = [curr.newImmHpvNonVax(1 : end , : , : , : , : , :); vaxResult{n}.newImmHpvNonVax(2 : end , : , : , : , : , :)];
vaxResult{n}.newScreen = [curr.newScreen(1 : end , :, :, :, :, :, :, :); vaxResult{n}.newScreen(2 : end , : , : , : , : , :, :, :)]; %[curr.newScreen(1 : end , : , : , : , : , : , : ); vaxResult{n}.newScreen(2 : end , : , : , : , : , : , :)];
vaxResult{n}.newHiv = [curr.newHiv(1 : end , : , : , : , : , : , :); vaxResult{n}.newHiv(2 : end , : , : , : , : , : , :)];
vaxResult{n}.hivDeaths = [curr.hivDeaths(1 : end , : , : , :); vaxResult{n}.hivDeaths(2 : end , : , : , :)];
% vaxResult{n}.artTreatTracker = [curr.artTreatTracker(1 : end , :  , : , : , : , :); vaxResult{n}.artTreatTracker(2 : end , : , : , : , : , :)];
vaxResult{n}.tVec = [curr.tVec(1 : end), vaxResult{n}.tVec(2 : end)];
% vaxResult{n}.ccSymp = [curr.ccSymp(1:end,:, :, :, :, :); vaxResult{n}.ccSymp(2:end,:, :, :, :, :)]; 
% vaxResult{n}.ccTreat = [curr.ccTreat(1:end, :, :, :, :, :); vaxResult{n}.ccTreat(2:end,:, :, :, :, :)]; 
vaxResult{n}.vaxdSchool = [curr.vaxdSchool(1:end, :); vaxResult{n}.vaxdSchool(2:end, :)]; % the only vax matrix in both historical and future sim 

% VACCINATIONS ********************************
    vaxTemplate = zeros(nTimepoints, 3); 
    vaxTemplate(1:end, 1) = vaxResult{n}.vaxdSchool; 
    vaxTemplate((size(curr.vaxdSchool,1)):end, 2) = vaxResult{n}.vaxdLmtd; % vaxdLmtd is only in futureSim, so add it to the point that futureSim starts
    vaxTemplate((size(curr.vaxdSchool,1)):end, 3) = vaxResult{n}.vaxdCU; % same as above
    vaxTotal = vaxTemplate(:, 1) + vaxTemplate(:, 2) + vaxTemplate(:, 3); % add up all the vaccine columns 

    vax(:, 17, 1, j) = vaxTemplate(:, 1); % input into the "all age" category 17; school based vax
    vax(:, 17, 2, j) = vaxTemplate(:, 3); % catchup vax 

% DEATHS (CC and all cause) **************************************
    ccDeath = zeros(nTimepoints, age); 
    for a = 1 : age 
        ccDeath_treat(:, a) = sum(sum(sum(vaxResult{n}.ccDeath_treat(:, :, a, :),2),3),4); 
        ccDeath_untreat(:, a) = sum(sum(sum(vaxResult{n}.ccDeath_untreat(:, :, a, :),2),3),4); 
    end 

    % combine all death data into 3D matrix
    deaths(:, 1:age, 1, j) = ccDeath_treat; % cc death stratified by age
    deaths(:, 1:age, 2, j) = ccDeath_untreat;
    deaths(:, 17, 3, j) = vaxResult{n}.deaths(:); % total all cause deaths not stratified by age

% CC HEALTH STATES (CC stage prevalence) ****************************

    for a = 1 : age 
        for x = 1 : endpoints 
            vaxInds1 = toInd(allcomb(1:disease, 1:viral, 1:hpvVaxStates, 6, x, 1:intervens, 2, a, 1:risk));
            vaxInds2 = toInd(allcomb(1:disease, 1:viral, 6, [1:5 7], x, 1:intervens, 2, a, 1:risk)); 
            % the purpose of [1:5 7] is to make sure there are no overlapping health states that you end up double counting
            vaxInds = [vaxInds1; vaxInds2]; 
            ccHealthState(1:end, a, x, j) = sum(vaxResult{n}.popVec(:, vaxInds), 2);
        end 
    end 

% HPV HEALTH STATES ************************************
% If you draw out a 7x7 matrix, you can map out what h and s values map out
% to an overarching hpv health state category. So for example, for h=5
% (CIN3), the overarching hpv health state is 5, so when h=5 and s is
% anything lower than that, or when s=5 and h is anything lower than that.
% when either s or h is 7 and the other is 5, that counts too. Note that I
% completely took out h or s = 6 and separate calculate the CC health
% states. I drew some pictures in by blue notebook that are helpful in
% developing this piece of code. 

    for a = 1 : age
        for h = 1 : hpvVaxStates
            if h < 7 
                vaxInds1 = toInd(allcomb(1:disease, 1:viral, [1:h 7], h, 1:endpoints, 1:intervens, 2, a, 1:risk)); 
                vaxInds2 = toInd(allcomb(1:disease, 1:viral, h, [1:(h-1) 7], 1:endpoints, 1:intervens, 2, a, 1:risk)); 
                vaxInds = [vaxInds1; vaxInds2]; 
                hpvHealthState(1:end, a, h, j) = sum(vaxResult{n}.popVec(:, vaxInds), 2);
            else % h=7 has different rules because only applicable is h=7 and s=7
                vaxInds = toInd(allcomb(1:disease, 1:viral, 7, 7, 1:endpoints, 1:intervens, 2, a, 1:risk)); 
                hpvHealthState(1:end, a, 7, j) = sum(vaxResult{n}.popVec(:, vaxInds), 2); % note 6. hpv immune will be at index 6 in hpvHealthState
            end 
        end
    end

% NEW CERVICAL CANCER CASES ******************************

    for a = 1 : age 
        newCC(:, a, j) = sum(sum(sum(vaxResult{n}.newCC(:, :, a, :),2),3),4);
    end 

% TOTAL NUMBER OF PEOPLE PER AGE GROUP ********************

    for a = 1 : age 
        vaxInds = toInd(allcomb(1:disease, 1:viral, 1:hpvVaxStates, 1:hpvNonVaxStates, 1:endpoints, 1:intervens, 2, a, 1:risk)); 
        totalPerAge(1:end, a, j) = sum(vaxResult{n}.popVec(:, vaxInds), 2); 
    end 

    % Adding index 17 to the age dimension for the total number of people of all ages 
    vaxInds = toInd(allcomb(1:disease, 1:viral, 1:hpvVaxStates, 1:hpvNonVaxStates, 1:endpoints, 1:intervens, 2, 1:age, 1:risk)); 
    totalPerAge(1:end, 17, j) = sum(vaxResult{n}.popVec(:, vaxInds), 2); 

% SCREENING ************************************************

% note that the results show 2 sceening age indices, but i believe results
% should only show up for one column of the matrix

colpoRetain = 0.72;
cinTreatRetain = 0.5;
ccTreatRetain = 0.4;

numScreen = zeros(nTimepoints, age, disease, hpvVaxStates, hpvNonVaxStates, 3); 
numColpo = numScreen; 
numCinTreat = numScreen; 
numCCTreat = zeros(nTimepoints, 3, age , 3, 2); % 2 is for treatment by screening or symptoms, and 3 is for treated, untreated, or hyst

    for aInd = 1 : length(screenAge) % modified to the ages that are screened  
%         a = sceScreenInd{aInd}; % the actual indices of newScreen  
        fullAgeInd = screenAge(aInd); % the indices (1:16) to actually put the results in
        a = screenAge(aInd); % indices within vaxResults.newScreen
%         fullAgeInd = sceScreenAges(aInd); % the indices (1:16) to actually put the results in
            for h = 1 : hpvVaxStates % Vaccine-type HPV precancer state
                for s = 1 : hpvNonVaxStates % Non-vaccine-type HPV precancer state
                    for x = 1 : 3 % Cervical cancer or hysterectomy status. we do not care to stratify by people with hyst, so endpoints-1. 
%                         for i = 3 : 4
                        % Apply selected screening/treatment algorithm
                            % if you're susceptible/immune to both HPV types or
                            % have had a hysterectomy 
                            if ( ((h<=3) || (h==7)) && ((s<=3) || (s==7)) && (x==1) ) || (x==10)

                                numScreen(: , fullAgeInd, h, s, x) = sum(sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , 1:disease , 1:viral , h , s , x , : , :),2),3),4),5),6),7),8);
                                numColpo(: , fullAgeInd, h, s, x) = 0.0; 
                                numCinTreat(: , fullAgeInd , h , s , x) = 0.0; 

                            % if you have CIN2+ of either HPV type
                            elseif ( (((h==4) || (h==5)) && ((s<=5) || (s==7))) || (((h<=5) || (h==7)) && ((s==4) || (s==5))) ) && (x==1) 
                                numScreen(: , fullAgeInd, h, s, x) = sum(sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , 1:disease , 1:viral , h , s , x , : , :),2),3),4),5),6),7),8);
                                numColpo(: , fullAgeInd , h , s , x) = numScreen(: , fullAgeInd , h , s , x) * screenAlgs{1}.testSens(2) * colpoRetain; 
                                numCinTreat(: , fullAgeInd , h , s , x) = numScreen(: , fullAgeInd , h , s , x) * screenAlgs{1}.testSens(2) * colpoRetain * cinTreatRetain; 

                            % if you have cervical cancer
                            elseif ( (x==1) && ((h==6) || (s==6)) ) || (x==2) || (x==3)
                                numScreen(: , fullAgeInd , h, s, x) = sum(sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , 1:disease , 1:viral , h , s , x , : , :),2),3),4),5),6),7),8);
                                numColpo(: , fullAgeInd , h , s , x) = numScreen(: , fullAgeInd , h , s , x) * screenAlgs{1}.testSens(2) * colpoRetain; 
                                numCinTreat(: , fullAgeInd , h , s , x) = 0.0; 

                            end
                  
                    end
                end
            end
    end

% DIAGNOSED CANCERS ************************************************

for a = 1:age
    for x = 1:3

        % val after x is 1 for treatment from screening or 2 for treatment from symptoms 
        % 2nd val is 1 for treated, 2 for untreated, and 3 for hyst
        % don't need to multiply by cc treat retain because LTFU is already accounted for within the model

        % edit for scenario 10: I condensed ccTreat and ccSymp for scenario 10 for future sim because of out of memory errors. 
        % below, i separately process historicalsim and futuresim results, and then add them together at the end
        % since the results have mismatched matrix sizes. 

        % curr - cc treatment from screening
        numCCTreat_curr(: , x , a , 1 , 1) = sum(sum(sum(sum(sum(curr.ccTreat(: , 1:disease , x , 1:4 , a , 1),2),3),4),5),6); % treated
        numCCTreat_curr(: , x , a , 2 , 1) = sum(sum(sum(sum(sum(curr.ccTreat(: , 1:disease , x , 1:4 , a , 2),2),3),4),5),6); % untreated
        numCCTreat_curr(: , x , a , 3 , 1) = sum(sum(sum(sum(sum(curr.ccTreat(: , 1:disease , x , 1:4 , a , 3),2),3),4),5),6); % hyst

        % curr - cc treatment from symptoms
        numCCTreat_curr(: , x , a , 1 , 2) = sum(sum(sum(sum(sum(curr.ccSymp(: , 1:disease , x , 1:4 , a , 1),2),3),4),5),6); % treated
        numCCTreat_curr(: , x , a , 2 , 2) = sum(sum(sum(sum(sum(curr.ccSymp(: , 1:disease , x , 1:4 , a , 2),2),3),4),5),6); % untreated
        numCCTreat_curr(: , x , a , 3 , 2) = sum(sum(sum(sum(sum(curr.ccSymp(: , 1:disease , x , 1:4 , a , 3),2),3),4),5),6); % hyst

        % fut - cc treatment from screening 
        numCCTreat_fut(: , x , a , 1 , 1) = sum(sum(sum(vaxResult{n}.ccTreat(: , x , a , 1),2),3),4); % treated
        numCCTreat_fut(: , x , a , 2 , 1) = sum(sum(sum(vaxResult{n}.ccTreat(: , x , a , 2),2),3),4); % untreated
        numCCTreat_fut(: , x , a , 3 , 1) = sum(sum(sum(vaxResult{n}.ccTreat(: , x , a , 3),2),3),4); % hyst

        % cc treatment from symptoms
        numCCTreat_fut(: , x , a , 1 , 2) = sum(sum(sum(vaxResult{n}.ccSymp(: , x , a , 1),2),3),4); % treated
        numCCTreat_fut(: , x , a , 2 , 2) = sum(sum(sum(vaxResult{n}.ccSymp(: , x , a , 2),2),3),4); % untreated
        numCCTreat_fut(: , x , a , 3 , 2) = sum(sum(sum(vaxResult{n}.ccSymp(: , x , a , 3),2),3),4); % hyst

    end 
end 

numCCTreat = [numCCTreat_curr; numCCTreat_fut(2:end, :, :, :, :)]; 


% SQUISHING SCREENING AND DIAGNOSED CANCERS ************************

    % Sum all dimensions so that we're only stratifying by time and age
    numScreenSquish = sum(sum(sum(numScreen(:, 1:age, 1:hpvVaxStates, 1:hpvNonVaxStates, 1:3), 3),4),5);
    numColpoSquish = sum(sum(sum(numColpo(:, 1:age, 1:hpvVaxStates, 1:hpvNonVaxStates, 1:3), 3),4),5);
    numCinTreatSquish = sum(sum(sum(numCinTreat(:, 1:age, 1:hpvVaxStates, 1:hpvNonVaxStates, 1:3), 3),4),5);

    % put screen/treat in a third dimension
    screenTreat(:, 1:age, 1, j) = numScreenSquish; 
    screenTreat(:, 1:age, 2, j) = numColpoSquish; 
    screenTreat(:, 1:age, 3, j) = numCinTreatSquish; 
    screenSympCCTreat(:, 1:3, 1:age , : , : , j) = numCCTreat; 

end % for loop end 

end % function end 