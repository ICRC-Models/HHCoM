function [deaths, screenTreat, hpvHealthState, ccHealthState, hivHealthState, totalPerAge, vax, nonDisabHealthState, scen0CinTx, newCC, screenSympCCTreat] = vaxCEA_multSims_SACEA_CH(vaxResultInd , sceNum , fileNameNums, fileInds, deaths, screenTreat, hpvHealthState, ccHealthState, hivHealthState,totalPerAge, vax, nonDisabHealthState, scen0CinTx, newCC, screenSympCCTreat)
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
    ccReghpvVaxInds_treat , ccDisthpvVaxInds_treat] = loadUp2(1 , 0 , [] , [] , [] , 1);

% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 1.5)

lastYear = 2123;  % ***SET ME***: last year of simulation (use 2122 for SA screening analysis)

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
baseFileName = ['22Apr20Ph2V112v57BaseVax_spCytoScreen_shortName_noVMMChpv_discontFxd_screenCovFxd_hivInt2017_SA-S' , sceNum , '_']; % ***SET ME***: name for simulation output file
% Looping length 
loopSegments = {0 , round(nRuns/2) , nRuns};
loopSegmentsLength = length(loopSegments);

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

if ((str2num(sceNum) == 0) || (str2num(sceNum) == 1)) % Scenarios 0 or 1
    % Screening paper cytology algorithm
    screenAlgs = spCyto;
    screenAlgs.genTypBool = 0; % whether or not method uses HPV genotyping, only looks for vaccine types, only get treated if you have the high risk types
    % proportion who return for treatment (susceptible/immune/infected/CIN1 , CIN2 , CIN3 , CC)
    screenAlgs.leepRetain = zeros(1,4);
    screenAlgs.cryoRetain = [0.0 , cryoRetain , cryoRetain , ccRetain];     
    screenAlgs.thrmlRetain = zeros(1,4);  
elseif ((str2num(sceNum) == 2) || (str2num(sceNum) == 3)) % Scenarios 2 or 3
    % Screening paper HPV DNA -and-treat algorithm
    screenAlgs = spHpvDna;
    screenAlgs.genTypBool = 0;
    % proportion who return for treatment (susceptible/immune/infected/CIN1 , CIN2 , CIN3 , CC)
    screenAlgs.leepRetain = [[eligLeep.*leepRetain] , ccRetain]; 
    screenAlgs.cryoRetain = zeros(1,4);
    screenAlgs.thrmlRetain = [[(1-eligLeep).*thrmlRetain] , ccRetain]; 
elseif ((str2num(sceNum) == 4) || (str2num(sceNum) == 5)) % Scenarios 4 or 5
    % Screening paper HPV DNA+genotyping -and-treat algorithm
    screenAlgs = spGentyp;
    screenAlgs.genTypBool = 1;
    % proportion who return for treatment (susceptible/immune/infected/CIN1 , CIN2 , CIN3 , CC)
    screenAlgs.leepRetain = [[eligLeep.*leepRetain] , ccRetain]; 
    screenAlgs.cryoRetain = zeros(1,4);
    screenAlgs.thrmlRetain = [[(1-eligLeep).*thrmlRetain] , ccRetain]; 
elseif ((str2num(sceNum) == 6) || (str2num(sceNum) == 7)) % Scenarios 6 or 7
    % Screening paper AVE -and-treat algorithm
    screenAlgs = spAve;
    screenAlgs.genTypBool = 0;
    % proportion who return for treatment (susceptible/immune/infected/CIN1 , CIN2 , CIN3 , CC)
    screenAlgs.leepRetain = [[eligLeep.*leepRetain] , ccRetain]; 
    screenAlgs.cryoRetain = zeros(1,4);
    screenAlgs.thrmlRetain = [[(1-eligLeep).*thrmlRetain] , ccRetain]; 
elseif ((str2num(sceNum) == 8) || (str2num(sceNum) == 9)) % Scenarios 8 or 9
    % Screening paper HPV DNA + AVE triage -and-treat algorithm
    screenAlgs = spHpvAve;
    screenAlgs.genTypBool = 0;
    % proportion who return for treatment (susceptible/immune/infected/CIN1 , CIN2 , CIN3 , CC)
    screenAlgs.leepRetain = [[eligLeep.*leepRetain] , ccRetain]; 
    screenAlgs.cryoRetain = zeros(1,4);
    screenAlgs.thrmlRetain = [[(1-eligLeep).*thrmlRetain] , ccRetain]; 
end

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
for k = 1 : loopSegmentsLength-1
    for j = loopSegments{k}+1 : loopSegments{k+1} % trying just a for loop for now

% k=1; % temporarily only running a single parameter to test 
% j=1; 

% Load results
pathModifier = [baseFileName , fileInds{j}];
nSims = size(dir([pwd , '/HHCoM_Results/' , pathModifier, '/' , '*.mat']) , 1);
curr = load([pwd , '/HHCoM_Results/toNow_22Apr20Ph2V112v57BaseVax_spCytoScreen_shortName_noVMMChpv_discontFxd_screenCovFxd_hivInt2017_' , fileInds{j}]); % ***SET ME***: name for historical run output file
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
vaxResult{n}.ccDeath_treat = [curr.ccDeath_treat(1 : end , : , : , :) ; vaxResult{n}.ccDeath_treat(2 : end , : , : , :)]; % consolidating historical CC death #s with future... etc.
vaxResult{n}.ccDeath_untreat = [curr.ccDeath_untreat(1 : end , : , : , :) ; vaxResult{n}.ccDeath_untreat(2 : end , : , : , :)]; 
vaxResult{n}.newCC = [curr.newCC(1 : end , : , : , :); vaxResult{n}.newCC(2 : end , : , : , :)];
vaxResult{n}.deaths = [curr.deaths(1 : end, 1); vaxResult{n}.deaths(2 : end, 1)];
vaxResult{n}.newHpvVax = [curr.newHpvVax(1 : end , : , : , : , : , :); vaxResult{n}.newHpvVax(2 : end , : , : , : , : , :)]; % infected with vaccine type HPV
vaxResult{n}.newImmHpvVax = [curr.newImmHpvVax(1 : end , : , : , : , : , :); vaxResult{n}.newImmHpvVax(2 : end , : , : , : , : , :)];
vaxResult{n}.newHpvNonVax = [curr.newHpvNonVax(1 : end , : , : , : , : , :); vaxResult{n}.newHpvNonVax(2 : end , : , : , : , : , :)];
vaxResult{n}.newImmHpvNonVax = [curr.newImmHpvNonVax(1 : end , : , : , : , : , :); vaxResult{n}.newImmHpvNonVax(2 : end , : , : , : , : , :)];
vaxResult{n}.newScreen = [vaxResult{n}.newScreen(1 : end , : , : , : , : , : , :)]; %[curr.newScreen(1 : end , : , : , : , : , : , : ); vaxResult{n}.newScreen(2 : end , : , : , : , : , : , :)];
vaxResult{n}.newHiv = [curr.newHiv(1 : end , : , : , : , : , : , :); vaxResult{n}.newHiv(2 : end , : , : , : , : , : , :)];
vaxResult{n}.hivDeaths = [curr.hivDeaths(1 : end , : , : , :); vaxResult{n}.hivDeaths(2 : end , : , : , :)];
vaxResult{n}.artTreatTracker = [curr.artTreatTracker(1 : end , :  , : , : , : , :); vaxResult{n}.artTreatTracker(2 : end , : , : , : , : , :)];
vaxResult{n}.tVec = [curr.tVec(1 : end), vaxResult{n}.tVec(2 : end)];
vaxResult{n}.vaxdSchool = [curr.vaxdSchool(:, :); vaxResult{n}.vaxdSchool(2:end, :)]; 
vaxResult{n}.ccSymp = [curr.ccSymp(1:end,:, :, :); vaxResult{n}.ccSymp(2:end,:, :, :)]; 
vaxResult{n}.ccTreat = [curr.ccTreat(1:end, :, :, :); vaxResult{n}.ccTreat(2:end,:, :, :)]; 

% VACCINATIONS ********************************
    vaxTemplate = zeros(nTimepoints, 3); 
    vaxTemplate(:, 1) = vaxResult{n}.vaxdSchool; 
    vaxTemplate((size(curr.vaxdSchool,1)):end, 2) = vaxResult{n}.vaxdLmtd; % vaxdLmtd is only in futureSim, so add it to the point that futureSim starts
    vaxTemplate((size(curr.vaxdSchool,1)):end, 3) = vaxResult{n}.vaxdCU; % same as above
    vaxTotal = vaxTemplate(:, 1) + vaxTemplate(:, 2) + vaxTemplate(:, 3); % add up all the vaccine columns 

    vax(:, 17, j) = vaxTotal; % input into the "all age" category 17

% DEATHS **************************************
    ccDeath = zeros(nTimepoints, age); 
    for a = 1 : age 
        ccDeath_treat(:, a) = sum(sum(sum(vaxResult{n}.ccDeath_treat(:, :, a, :),2),3),4); 
        ccDeath_untreat(:, a) = sum(sum(sum(vaxResult{n}.ccDeath_untreat(:, :, a, :),2),3),4); 
    end 

    hivDeath = zeros(nTimepoints, age); 
    for a = 1 : age 
        hivDeath(:, a) = sum(sum(sum(vaxResult{n}.hivDeaths(:, :, 2, 1),2),3),4); 
    end 

    % combine all death data into 3D matrix
    deaths(:, 1:age, 1, j) = ccDeath_treat; % cc death stratified by age
    deaths(:, 1:age, 2, j) = ccDeath_untreat;
    deaths(:, 1:age, 3, j) = hivDeath; 
    deaths(:, 17, 4, j) = vaxResult{n}.deaths(:); 

% CC HEALTH STATES ****************************

    for a = 1 : age 
	    for x = 1 : endpoints 
		    vaxInds1 = toInd(allcomb(1:disease, 1:viral, 1:hpvVaxStates, 6, x, 1:intervens, 2, a, 1:risk));
            vaxInds2 = toInd(allcomb(1:disease, 1:viral, 6, [1:5 7], x, 1:intervens, 2, a, 1:risk)); 
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

% NON DISABILITY HEALTH STATES ***************************

    nonDisabVector = [1 2 3 4 5 7]; % indices for the h and s comparments for the non-disability health states (everything except for cervical cancer or hysterectomy)

    for a = 1 : age 
        % pull popVec indices for all non-disability health states
        % h and s are based on nonDisabVector
        % disease is indices 1 and 2 (HIV negative)
        nonDisabInds = toInd(allcomb(1:2, 1:viral, nonDisabVector, nonDisabVector, 1:endpoints, 1:intervens, 2, a, 1:risk)); 
        nonDisabHealthState(1:end, a, j) = sum(vaxResult{n}.popVec(:, nonDisabInds), 2);
    end 

% HIV HEALTH STATES ************************************

    for a = 1 : age
        for dInd = 1 : length(diseaseVec_vax)
            d = diseaseVec_vax{dInd}; 

            vaxInds = toInd(allcomb(d, 1:viral, 1:hpvVaxStates, 1:hpvNonVaxStates, 1:endpoints, 1:intervens, 2, a, 1:risk)); 
            hivHealthState(1:end, a, dInd, j) = sum(vaxResult{n}.popVec(:, vaxInds), 2); 
        end
    end 

% NEW CERVICAL CANCER CASES ******************************

    % vaxResult{n}.newCC % (time, disease, age, hpvTypeGroups (2))

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

numScreen = zeros(nTimepointsScreen, disease, hpvVaxStates, hpvNonVaxStates, 3, age); 
numLEEP = zeros(nTimepointsScreen, disease, hpvVaxStates, hpvNonVaxStates, 3, age); 
numCryo = zeros(nTimepointsScreen, disease, hpvVaxStates, hpvNonVaxStates, 3, age); 
numThrml = zeros(nTimepointsScreen, disease, hpvVaxStates, hpvNonVaxStates, 3, age); 
numHyst = zeros(nTimepointsScreen, disease, hpvVaxStates, hpvNonVaxStates, 3, age); 

    for aInd = 1 : length(sceScreenInd) % modified to the ages that are screened  
%         a = sceScreenInd{aInd}; % the actual indices of newScreen  
        fullAgeInd = sceScreenIndAge(aInd); % the indices (1:16) to actually put the results in
        a = sceScreenInd{aInd}; % indices within vaxResults.newScreen
%         fullAgeInd = sceScreenAges(aInd); % the indices (1:16) to actually put the results in
        for dInd = 1 : disease % HIV disease states in template
            d = dInd; % originally d was separate from dInd due to the combination of d=2 and d=1 as HIV neg. now we don't care about stratifying. 
            for h = 1 : hpvVaxStates % Vaccine-type HPV precancer state
                for s = 1 : hpvNonVaxStates % Non-vaccine-type HPV precancer state
                    for x = 1 : 3 % Cervical cancer or hysterectomy status. we do not care to stratify by people with hyst, so endpoints-1. 
                        % Apply selected screening/treatment algorithm
                            % if you're susceptible/immune to both HPV types or
                            % have had a hysterectomy 
                            if [( ((h==1) || (h==7)) && ((s==1) || (s==7)) && (x==1) )] || (x==4) || [(screenAlgs.genTypBool && ((h==1) || (h==7)) && (((s>=2) && (s<=5)) || ((s==6) && (x<=3))))]
                                numScreen(: , dInd , h, s, x, fullAgeInd) = sum(sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , d , 1:viral, h , s , x , a , :),2),3),4),5),6),7),8);
                                numLEEP(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * 0.0;
                                numCryo(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * 0.0;
                                numThrml(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * 0.0;
                                numHyst(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * 0.0;

                            % if you're infected with either HPV type
                            elseif [( ((h==2) && ((s<=2) || (s==7))) || (((h<=2) || (h==7)) && (s==2)) ) && (x==1)] && [(~screenAlgs.genTypBool) || (screenAlgs.genTypBool && (h==2))]
                                numScreen(: , dInd , h, s, x, fullAgeInd) = sum(sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , d , 1:viral, h , s , x , a , :),2),3),4),5),6),7),8);
                                numLEEP(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * screenAlgs.testSens(d,2) * screenAlgs.colpoRetain * screenAlgs.leepRetain(1);
                                numCryo(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * screenAlgs.testSens(d,2) * screenAlgs.colpoRetain * screenAlgs.cryoRetain(1);
                                numThrml(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * screenAlgs.testSens(d,2) * screenAlgs.colpoRetain * screenAlgs.thrmlRetain(1);
                                numHyst(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * 0.0;

                            % if you have CIN1 of either HPV type
                            elseif [( ((h==3) && ((s<=3) || (s==7))) || (((h<=3) || (h==7)) && (s==3)) ) && (x==1)] && [(~screenAlgs.genTypBool) || (screenAlgs.genTypBool && (h==3))]
                                numScreen(: , dInd , h, s, x, fullAgeInd) = sum(sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , d , 1:viral, h , s , x , a , :),2),3),4),5),6),7),8);
                                numLEEP(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * screenAlgs.testSens(d,2) * screenAlgs.colpoRetain * screenAlgs.leepRetain(1);
                                numCryo(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * screenAlgs.testSens(d,2) * screenAlgs.colpoRetain * screenAlgs.cryoRetain(1);
                                numThrml(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * screenAlgs.testSens(d,2) * screenAlgs.colpoRetain * screenAlgs.thrmlRetain(1);
                                numHyst(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * 0.0;

                            % if you have CIN2 of either HPV type
                            elseif [( ((h==4) && ((s<=4) || (s==7))) || (((h<=4) || (h==7)) && (s==4)) ) && (x==1)] && [(~screenAlgs.genTypBool) || (screenAlgs.genTypBool && (h==4))]
                                numScreen(: , dInd , h, s, x, fullAgeInd) = sum(sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , d , 1:viral, h , s , x , a , :),2),3),4),5),6),7),8);
                                numLEEP(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * screenAlgs.testSens(d,3) * screenAlgs.colpoRetain * screenAlgs.leepRetain(2);
                                numCryo(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * screenAlgs.testSens(d,3) * screenAlgs.colpoRetain * screenAlgs.cryoRetain(2);
                                numThrml(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * screenAlgs.testSens(d,3) * screenAlgs.colpoRetain * screenAlgs.thrmlRetain(2);
                                numHyst(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * 0.0;
                            % if you have CIN3 of either HPV type
                            elseif [( ((h==5) && ((s<=5) || (s==7))) || (((h<=5) || (h==7)) && (s==5)) ) && (x==1)] && [(~screenAlgs.genTypBool) || (screenAlgs.genTypBool && (h==5))]
                                numScreen(: , dInd , h, s, x, fullAgeInd) = sum(sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , d , 1:viral, h , s , x , a , :),2),3),4),5),6),7),8);
                                numLEEP(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * screenAlgs.testSens(d,4) * screenAlgs.colpoRetain * screenAlgs.leepRetain(3);
                                numCryo(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * screenAlgs.testSens(d,4) * screenAlgs.colpoRetain * screenAlgs.cryoRetain(3);
                                numThrml(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * screenAlgs.testSens(d,4) * screenAlgs.colpoRetain * screenAlgs.thrmlRetain(3);
                                numHyst(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * 0.0;
                            % if you have cervical cancer
                            elseif [( (x==1) && ((h==6) || (s==6)) ) || (x==2) || (x==3)] && [(~screenAlgs.genTypBool) || (screenAlgs.genTypBool && (h==6) && (x<=3))]
                                numScreen(: , dInd , h, s, x, fullAgeInd) = sum(sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , d , 1:viral, h , s , x , a , :),2),3),4),5),6),7),8);
                                numLEEP(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * 0.0;
                                numCryo(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * 0.0;
                                numThrml(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * 0.0;
                                numHyst(: , dInd , h, s, x, fullAgeInd) = numScreen(: , dInd , h, s, x, fullAgeInd) * screenAlgs.testSens(d,4) * screenAlgs.colpoRetain * screenAlgs.treatRetain(4);
                            end
                    end
                end
            end
        end
    end

    % Sum all dimensions so that we're only stratifying by time and age
    numScreenSquish = sum(sum(sum(sum(numScreen(:, 1:disease, 1:hpvVaxStates, 1:hpvNonVaxStates, 1:3, 1:age), 2),3),4),5);
    numLEEPSquish = sum(sum(sum(sum(numLEEP(:, 1:disease, 1:hpvVaxStates, 1:hpvNonVaxStates, 1:3, 1:age), 2),3),4),5);
    numCryoSquish = sum(sum(sum(sum(numCryo(:, 1:disease, 1:hpvVaxStates, 1:hpvNonVaxStates, 1:3, 1:age), 2),3),4),5);
    numThrmlSquish = sum(sum(sum(sum(numThrml(:, 1:disease, 1:hpvVaxStates, 1:hpvNonVaxStates, 1:3, 1:age), 2),3),4),5);
    numHystSquish = sum(sum(sum(sum(numHyst(:, 1:disease, 1:hpvVaxStates, 1:hpvNonVaxStates, 1:3, 1:age), 2),3),4),5);

    % put screen/treat in a third dimension
    screenTreat(:, 1:age, 1, j) = numScreenSquish; 
    screenTreat(:, 1:age, 2, j) = numLEEPSquish; 
    screenTreat(:, 1:age, 3, j) = numCryoSquish; 
    screenTreat(:, 1:age, 4, j) = numThrmlSquish; 
    screenTreat(:, 1:age, 5, j) = numHystSquish; 

    % DIAGNOSED CANCERS ************************************************

    for a = 1:age
        for x = 1:3

            % val after x is 1 for treatment from screening or 2 for treatment from symptoms 
            % 2nd val is 1 for treated, 2 for untreated, and 3 for hyst
            % don't need to multiply by cc treat retain because LTFU is already accounted for within the model

            % cc treatment from screening 
            numCCTreat(: , x , a , 1 , 1) = sum(sum(sum(vaxResult{n}.ccTreat(: , x , a , 1),2),3),4); % treated
            numCCTreat(: , x , a , 2 , 1) = sum(sum(sum(vaxResult{n}.ccTreat(: , x , a , 2),2),3),4); % untreated
            numCCTreat(: , x , a , 3 , 1) = sum(sum(sum(vaxResult{n}.ccTreat(: , x , a , 3),2),3),4); % hyst

            % cc treatment from symptoms
            numCCTreat(: , x , a , 1 , 2) = sum(sum(sum(vaxResult{n}.ccSymp(: , x , a , 1),2),3),4); % treated
            numCCTreat(: , x , a , 2 , 2) = sum(sum(sum(vaxResult{n}.ccSymp(: , x , a , 2),2),3),4); % untreated
            numCCTreat(: , x , a , 3 , 2) = sum(sum(sum(vaxResult{n}.ccSymp(: , x , a , 3),2),3),4); % hyst

        end 
    end 

    % SQUISHING SCREENING AND DIAGNOSED CANCERS ************************

    screenSympCCTreat(:, 1:3, 1:age , : , : , j) = numCCTreat; 

    % Separate therapy count stratified by CC state for a one-off analysis
    % Jacinda wants to do in her CEA only for scenario 0
    % I used the same method of isolating h based on the highest value
    % between the vax and nonvax hpv health state
    for h = 1 : hpvVaxStates
        if h < 7 
            numLEEP_1 = sum(sum(sum(sum(sum(numLEEP(:, 1:disease, [1:h 7], h, 1:3, 1:age),2),3),4),5),6);
            numLEEP_2 = sum(sum(sum(sum(sum(numLEEP(:, 1:disease, h, [1:(h-1) 7], 1:3, 1:age),2),3),4),5),6);
            numLEEP_scen0 = numLEEP_1 + numLEEP_2; % hopefully this should result in a vector ??
            
            numCryo_1 = sum(sum(sum(sum(sum(numCryo(:, 1:disease, [1:h 7], h, 1:3, 1:age),2),3),4),5),6);
            numCryo_2 = sum(sum(sum(sum(sum(numCryo(:, 1:disease, h, [1:(h-1) 7], 1:3, 1:age),2),3),4),5),6);
            numCryo_scen0 = numCryo_1 + numCryo_2; % hopefully this should result in a vector ??

            scen0CinTx(:, h, 1, j) = numLEEP_scen0; 
            scen0CinTx(:, h, 2, j) = numCryo_scen0; 
        else % h=7 has different rules because only applicable is h=7 and s=7
            numLEEP_scen0 = sum(sum(sum(sum(sum(numLEEP(:, 1:disease, 7, 7, 1:3, 1:age),2),3),4),5),6); 
            numCryo_scen0 = sum(sum(sum(sum(sum(numCryo(:, 1:disease, 7, 7, 1:3, 1:age),2),3),4),5),6); 

            scen0CinTx(:, h, 1, j) = numLEEP_scen0; 
            scen0CinTx(:, h, 2, j) = numCryo_scen0; 
        end 
    end


    end
end % for loop end 

end % function end 